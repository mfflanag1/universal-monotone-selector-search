#!/usr/bin/env python3
"""Search for K2,3 terminal envelopes by local-game cutting planes."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
sys.path.append(str(PROJECT_SRC))

from n5_facet_search import load_facets
from n5_k23_batch import permute_mask
from n5_middle_split_identity_certificate import coefficient_vector
from n5_mixed_divergence_premium import choquet_representation
from n5_mixed_terminal_branch_search import extreme_representations
from n5_perspective_envelope_certificate import (
    embedded_topology_generators,
    protected_pairs_from_report,
)


F = Fraction


def sparse_matrix(rows: list[dict[int, float]], columns: int):
    row_indices = []
    column_indices = []
    values = []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                row_indices.append(row_index)
                column_indices.append(column)
                values.append(value)
    return coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(rows), columns),
    ).tocsr()


def permute_vector(
    vector: tuple[F, ...], permutation: tuple[int, ...]
) -> tuple[F, ...]:
    result = [F(0)] * len(vector)
    for player, value in enumerate(vector):
        result[permutation[player]] = value
    return tuple(result)


def permute_coefficients(
    coefficients: tuple[F, ...], permutation: tuple[int, ...]
) -> tuple[F, ...]:
    result = [F(0)] * len(coefficients)
    for coalition, value in enumerate(coefficients):
        result[permute_mask(coalition, permutation)] = value
    return tuple(result)


def seed_envelopes(
    reference_path: Path,
    divergences: dict[int, tuple[F, ...]],
) -> dict[int, tuple[F, ...]]:
    reference = json.loads(reference_path.read_text())
    reference_divergences = {
        int(node): tuple(F(value) for value in values)
        for node, values in reference["terminal_divergences"].items()
    }
    reference_envelopes = {
        int(node): tuple(F(value) for value in values)
        for node, values in reference["envelopes"].items()
    }
    result = {}
    for node, target in divergences.items():
        donors = range(2) if node < 2 else range(2, 5)
        chosen = None
        for donor in donors:
            for permutation in itertools.permutations(range(len(target))):
                if permute_vector(
                    reference_divergences[donor], permutation
                ) == target:
                    chosen = permute_coefficients(
                        reference_envelopes[donor], permutation
                    )
                    break
            if chosen is not None:
                break
        if chosen is None:
            # A reference library need not contain every divergence orbit.
            # The seed affects only the L1 tie-breaking objective, not the
            # topology or local-validity constraints, so use the canonical
            # Choquet branch as a deterministic fallback for a new type.
            chosen = coefficient_vector(
                choquet_representation(target, len(target)), len(target)
            )
        result[node] = chosen
    return result


def normalized_game_cone(n: int):
    grand = (1 << n) - 1
    variable_count = grand + 1
    facets, _ = load_facets()
    inequalities: list[dict[int, float]] = []
    rhs = []
    for facet in facets:
        inequalities.append(
            {
                coalition: -float(coefficient)
                for coalition, coefficient in enumerate(facet)
                if coefficient
            }
        )
        rhs.append(0.0)
    for coalition in range(variable_count):
        for player in range(n):
            if coalition >> player & 1:
                continue
            successor = coalition | (1 << player)
            inequalities.append({coalition: 1.0, successor: -1.0})
            rhs.append(0.0)
    equalities = np.zeros((2, variable_count))
    equalities[0, 0] = 1.0
    equalities[1, grand] = 1.0
    return (
        sparse_matrix(inequalities, variable_count),
        np.asarray(rhs),
        equalities,
        np.asarray([0.0, 1.0]),
    )


def solve_master(
    seed: dict[int, tuple[F, ...]],
    divergences: dict[int, tuple[F, ...]],
    topology_generators: list[tuple[F, ...]],
    cuts: list[tuple[int, int, tuple[float, ...], float]],
) -> tuple[dict[int, np.ndarray], float]:
    node_ids = sorted(divergences)
    n = len(divergences[node_ids[0]])
    coalition_count = 1 << n
    envelope_count = len(node_ids) * coalition_count
    topology_offset = envelope_count
    deviation_offset = topology_offset + len(topology_generators)
    variable_count = deviation_offset + envelope_count

    equalities: list[dict[int, float]] = []
    equality_rhs = []
    for coordinate in range(envelope_count):
        row = {coordinate: 1.0}
        for generator, topology_column in zip(
            topology_generators,
            range(topology_offset, deviation_offset),
            strict=True,
        ):
            if generator[coordinate]:
                row[topology_column] = float(generator[coordinate])
        equalities.append(row)
        equality_rhs.append(0.0)
    for node_offset, node in enumerate(node_ids):
        for player, target in enumerate(divergences[node]):
            equalities.append(
                {
                    node_offset * coalition_count + coalition: 1.0
                    for coalition in range(coalition_count)
                    if coalition >> player & 1
                }
            )
            equality_rhs.append(float(target))

    inequalities: list[dict[int, float]] = []
    inequality_rhs = []
    for node_offset, node in enumerate(node_ids):
        for coalition in range(coalition_count):
            coordinate = node_offset * coalition_count + coalition
            deviation = deviation_offset + coordinate
            expected = float(seed[node][coalition])
            inequalities.append({coordinate: 1.0, deviation: -1.0})
            inequality_rhs.append(expected)
            inequalities.append({coordinate: -1.0, deviation: -1.0})
            inequality_rhs.append(-expected)
    for node, _, game, branch_value in cuts:
        node_offset = node_ids.index(node) * coalition_count
        inequalities.append(
            {
                node_offset + coalition: -value
                for coalition, value in enumerate(game)
                if value
            }
        )
        inequality_rhs.append(-branch_value)

    objective = np.zeros(variable_count)
    objective[deviation_offset:] = 1.0
    bounds = (
        [(-20.0, 20.0)] * envelope_count
        + [(0.0, None)] * len(topology_generators)
        + [(0.0, None)] * envelope_count
    )
    solved = linprog(
        objective,
        A_ub=sparse_matrix(inequalities, variable_count),
        b_ub=np.asarray(inequality_rhs),
        A_eq=sparse_matrix(equalities, variable_count),
        b_eq=np.asarray(equality_rhs),
        bounds=bounds,
        method="highs-ds",
        options={
            "dual_feasibility_tolerance": 1e-10,
            "primal_feasibility_tolerance": 1e-10,
        },
    )
    if not solved.success:
        raise RuntimeError(f"envelope master failed: {solved.message}")
    primary_solution = solved
    primary_distance = float(solved.fun)
    inequalities.append(
        {
            deviation_offset + coordinate: 1.0
            for coordinate in range(envelope_count)
        }
    )
    inequality_rhs.append(primary_distance + 1e-9)
    secondary_objective = np.zeros(variable_count)
    for coordinate in range(envelope_count):
        secondary_objective[coordinate] = (
            ((coordinate * 37) % 101) - 50
        ) / 50.0
    for topology_column in range(topology_offset, deviation_offset):
        secondary_objective[topology_column] = (
            ((topology_column * 17) % 29) + 1
        ) * 1e-8
    selected = linprog(
        secondary_objective,
        A_ub=sparse_matrix(inequalities, variable_count),
        b_ub=np.asarray(inequality_rhs),
        A_eq=sparse_matrix(equalities, variable_count),
        b_eq=np.asarray(equality_rhs),
        bounds=bounds,
        method="highs-ds",
        options={
            "dual_feasibility_tolerance": 1e-10,
            "primal_feasibility_tolerance": 1e-10,
        },
    )
    if selected.success:
        solved = selected
    else:
        # The secondary objective only selects a stable point on the
        # primary-optimal face.  HiGHS can occasionally report Unknown on
        # this deliberately thin face even though the primary optimum is a
        # valid master solution.  Falling back preserves every mathematical
        # constraint and lets separation provide the next exact cut.
        solved = primary_solution
    envelopes = {
        node: solved.x[
            node_offset * coalition_count : (node_offset + 1) * coalition_count
        ]
        for node_offset, node in enumerate(node_ids)
    }
    return envelopes, primary_distance


def separating_cuts(
    envelopes: dict[int, np.ndarray],
    representations: dict[int, list[tuple[F, tuple[tuple[int, F], ...]]]],
    n: int,
    game_inequalities,
    game_rhs,
    game_equalities,
    game_equality_rhs,
    tolerance: float,
    per_node: int,
) -> tuple[list[tuple[int, int, tuple[float, ...], float]], float]:
    violations = []
    for node, branches in representations.items():
        node_violations = []
        envelope = envelopes[node]
        for branch, representation in enumerate(branches):
            coefficients = np.asarray(
                [float(value) for value in coefficient_vector(representation, n)]
            )
            difference = coefficients - envelope
            solved = linprog(
                -difference,
                A_ub=game_inequalities,
                b_ub=game_rhs,
                A_eq=game_equalities,
                b_eq=game_equality_rhs,
                bounds=[(None, None)] * len(difference),
                method="highs-ds",
                options={
                    "dual_feasibility_tolerance": 1e-10,
                    "primal_feasibility_tolerance": 1e-10,
                },
            )
            if not solved.success:
                raise RuntimeError("local envelope separation LP failed")
            violation = -float(solved.fun)
            if violation > tolerance:
                branch_value = float(coefficients @ solved.x)
                node_violations.append(
                    (
                        violation,
                        node,
                        branch,
                        tuple(float(value) for value in solved.x),
                        branch_value,
                    )
                )
        node_violations.sort(reverse=True)
        violations.extend(node_violations[:per_node])
    violations.sort(reverse=True)
    return (
        [
            (node, branch, game, branch_value)
            for _, node, branch, game, branch_value in violations
        ],
        max((row[0] for row in violations), default=0.0),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("--seed-reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--cuts-per-node", type=int, default=20)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    terminals = report["terminals"]
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    divergences = {
        int(terminal["node"]): tuple(
            F(value) for value in terminal["scaled_divergence"]
        )
        for terminal in terminals
    }
    n = len(next(iter(divergences.values())))
    representations = {
        node: extreme_representations(list(divergence), n)
        for node, divergence in divergences.items()
    }
    protected_pairs = protected_pairs_from_report(report)
    topology_generators, _ = embedded_topology_generators(
        node_ids, n, protected_pairs
    )
    seed = seed_envelopes(args.seed_reference, divergences)
    game_cone = normalized_game_cone(n)
    cuts: list[tuple[int, int, tuple[float, ...], float]] = []
    seen_cuts = set()
    envelopes = None
    maximum_violation = None
    for iteration in range(args.iterations):
        envelopes, master_distance = solve_master(
            seed, divergences, topology_generators, cuts
        )
        new_cuts, maximum_violation = separating_cuts(
            envelopes,
            representations,
            n,
            *game_cone,
            args.tolerance,
            args.cuts_per_node,
        )
        added = 0
        for node, branch, game, branch_value in new_cuts:
            key = (
                node,
                tuple(round(value, 10) for value in game),
            )
            if key in seen_cuts:
                continue
            seen_cuts.add(key)
            cuts.append((node, branch, game, branch_value))
            added += 1
        print(
            json.dumps(
                {
                    "iteration": iteration,
                    "master_l1_distance": master_distance,
                    "cut_count": len(cuts),
                    "new_cut_count": added,
                    "maximum_local_violation": maximum_violation,
                }
            ),
            flush=True,
        )
        if maximum_violation <= args.tolerance:
            break
        if added == 0:
            raise RuntimeError("cut generation stalled above tolerance")
    if envelopes is None or maximum_violation is None:
        raise RuntimeError("CEGAR did not run")
    payload = {
        "status": (
            "n5_k23_envelope_cegar_converged"
            if maximum_violation <= args.tolerance
            else "n5_k23_envelope_cegar_incomplete"
        ),
        "source": str(args.terminal_report),
        "seed_reference": str(args.seed_reference),
        "iterations": iteration + 1,
        "cut_count": len(cuts),
        "maximum_local_violation": maximum_violation,
        "cuts": [
            {
                "node": node,
                "branch": branch,
                "game": [format(value, ".17g") for value in game],
                "branch_value": format(branch_value, ".17g"),
            }
            for node, branch, game, branch_value in cuts
        ],
        "active_cuts": [
            {
                "node": node,
                "branch": branch,
                "game": [format(value, ".17g") for value in game],
                "branch_value": format(branch_value, ".17g"),
                "slack": format(
                    float(envelopes[node] @ np.asarray(game) - branch_value),
                    ".17g",
                ),
            }
            for node, branch, game, branch_value in cuts
            if float(envelopes[node] @ np.asarray(game) - branch_value)
            <= 1e-7
        ],
        "envelopes": {
            str(node): [format(float(value), ".17g") for value in envelope]
            for node, envelope in envelopes.items()
        },
        "perspective_link_support": {
            "inequalities": [],
            "equalities": [
                {
                    "row": ["perspective_link", node, coalition],
                    "weight": format(float(value), ".17g"),
                }
                for node, envelope in envelopes.items()
                for coalition, value in enumerate(envelope)
                if abs(float(value)) > 1e-14
            ],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0 if maximum_violation <= args.tolerance else 1


if __name__ == "__main__":
    raise SystemExit(main())
