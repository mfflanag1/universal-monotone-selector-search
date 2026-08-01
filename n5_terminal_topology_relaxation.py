#!/usr/bin/env python3
"""Optimize an indicator-dual objective using only terminal path implications."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csr_matrix

from n5_facet_search import load_facets


F = Fraction


def sparse_matrix(
    rows: list[dict[int, float]], column_count: int
) -> csr_matrix:
    row_indices = []
    columns = []
    values = []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                row_indices.append(row_index)
                columns.append(column)
                values.append(value)
    return coo_matrix(
        (values, (row_indices, columns)),
        shape=(len(rows), column_count),
    ).tocsr()


def reconstruct_exact_dual(
    result: Any,
    inequalities: list[dict[int, float]],
    inequality_rhs: list[float],
    inequality_metadata: list[tuple[Any, ...]],
    equalities: list[dict[int, float]],
    equality_rhs: list[float],
    equality_metadata: list[tuple[Any, ...]],
    objective: list[F],
    objective_constant: F | None,
    active_tolerance: float = 1e-9,
) -> dict[str, Any] | None:
    active_inequalities = [
        index
        for index, value in enumerate(result.ineqlin.marginals)
        if abs(float(value)) > active_tolerance
    ]
    active_equalities = [
        index
        for index, value in enumerate(result.eqlin.marginals)
        if abs(float(value)) > active_tolerance
    ]

    def validate(
        inequality_weights: dict[int, F],
        equality_weights: dict[int, F],
        max_denominator: int,
        reconstruction_method: str,
    ) -> dict[str, Any] | None:
        inequality_weights = {
            index: weight
            for index, weight in inequality_weights.items()
            if weight
        }
        equality_weights = {
            index: weight for index, weight in equality_weights.items() if weight
        }
        if any(weight > 0 for weight in inequality_weights.values()):
            return None
        stationarity = [F(0)] * len(objective)
        for index, weight in inequality_weights.items():
            for column, coefficient in inequalities[index].items():
                stationarity[column] += weight * F(coefficient)
        for index, weight in equality_weights.items():
            for column, coefficient in equalities[index].items():
                stationarity[column] += weight * F(coefficient)
        if stationarity != objective:
            return None
        dual_objective = sum(
            (
                weight * F(inequality_rhs[index])
                for index, weight in inequality_weights.items()
            ),
            F(0),
        ) + sum(
            (
                weight * F(equality_rhs[index])
                for index, weight in equality_weights.items()
            ),
            F(0),
        )
        if (
            objective_constant is not None
            and dual_objective != objective_constant
        ):
            return None
        return {
            "reconstruction_max_denominator": max_denominator,
            "reconstruction_method": reconstruction_method,
            "variable_objective": str(dual_objective),
            "lower_expectation_objective_with_constant": (
                None
                if objective_constant is None
                else str(-dual_objective + objective_constant)
            ),
            "active_inequalities": [
                {
                    "row": list(inequality_metadata[index]),
                    "weight": str(weight),
                }
                for index, weight in inequality_weights.items()
            ],
            "active_equalities": [
                {
                    "row": list(equality_metadata[index]),
                    "weight": str(weight),
                }
                for index, weight in equality_weights.items()
            ],
        }

    denominators = (
        1_000,
        10_000,
        100_000,
        1_000_000,
        10_000_000,
        100_000_000,
        1_000_000_000,
        10_000_000_000,
        1_000_000_000_000,
    )
    for max_denominator in denominators:
        certificate = validate(
            {
                index: F(
                    float(result.ineqlin.marginals[index])
                ).limit_denominator(max_denominator)
                for index in active_inequalities
            },
            {
                index: F(
                    float(result.eqlin.marginals[index])
                ).limit_denominator(max_denominator)
                for index in active_equalities
            },
            max_denominator,
            "direct_marginal_reconstruction",
        )
        if certificate is not None:
            return certificate

    support = [
        ("inequality", index) for index in active_inequalities
    ] + [("equality", index) for index in active_equalities]
    row_indices = []
    columns = []
    values = []
    for support_column, (row_type, index) in enumerate(support):
        if row_type == "inequality":
            row = inequalities[index]
            rhs = inequality_rhs[index]
        else:
            row = equalities[index]
            rhs = equality_rhs[index]
        for variable, coefficient in row.items():
            row_indices.append(variable)
            columns.append(support_column)
            values.append(coefficient)
        if rhs and objective_constant is not None:
            row_indices.append(len(objective))
            columns.append(support_column)
            values.append(rhs)
    repair_row_count = len(objective) + int(objective_constant is not None)
    repair_matrix = coo_matrix(
        (values, (row_indices, columns)),
        shape=(repair_row_count, len(support)),
    ).tocsr()
    repair_target = list(objective)
    if objective_constant is not None:
        repair_target.append(objective_constant)
    repair_rhs = np.asarray([float(value) for value in repair_target])
    repair_bounds = [
        (None, 0.0) if row_type == "inequality" else (None, None)
        for row_type, _ in support
    ]

    def support_residual(weights: list[F]) -> list[F]:
        lhs = [F(0)] * repair_row_count
        for support_column, (row_type, index) in enumerate(support):
            weight = weights[support_column]
            if not weight:
                continue
            if row_type == "inequality":
                row = inequalities[index]
                rhs = inequality_rhs[index]
            else:
                row = equalities[index]
                rhs = equality_rhs[index]
            for variable, coefficient in row.items():
                lhs[variable] += weight * F(coefficient)
            if objective_constant is not None:
                lhs[-1] += weight * F(rhs)
        target = repair_target
        return [
            expected - actual
            for expected, actual in zip(target, lhs, strict=True)
        ]

    current = [
        F(
            float(
                result.ineqlin.marginals[index]
                if row_type == "inequality"
                else result.eqlin.marginals[index]
            )
        ).limit_denominator(1_000_000)
        for row_type, index in support
    ]
    for refinement in range(12):
        residual = support_residual(current)
        scale = max((abs(value) for value in residual), default=F(0))
        if scale == 0:
            break
        correction = linprog(
            np.zeros(len(support)),
            A_eq=repair_matrix,
            b_eq=np.asarray([float(value / scale) for value in residual]),
            bounds=[(None, None)] * len(support),
            method="highs-ds",
        )
        if not correction.success:
            break
        best_candidate = None
        best_error = None
        for max_denominator in denominators:
            candidate = [
                value
                + scale
                * F(float(delta)).limit_denominator(max_denominator)
                for value, delta in zip(
                    current, correction.x, strict=True
                )
            ]
            inequality_weights = {
                index: candidate[support_column]
                for support_column, (row_type, index) in enumerate(support)
                if row_type == "inequality"
            }
            equality_weights = {
                index: candidate[support_column]
                for support_column, (row_type, index) in enumerate(support)
                if row_type == "equality"
            }
            certificate = validate(
                inequality_weights,
                equality_weights,
                max_denominator,
                f"rational_residual_refinement_{refinement}",
            )
            if certificate is not None:
                return certificate
            candidate_residual = support_residual(candidate)
            error = max(
                (float(abs(value)) for value in candidate_residual),
                default=0.0,
            )
            if best_error is None or error < best_error:
                best_error = error
                best_candidate = candidate
        assert best_candidate is not None
        current = best_candidate

    rng = np.random.default_rng(20260731)
    repair_objectives = [np.zeros(len(support))]
    repair_objectives.extend(
        rng.integers(-5, 6, size=len(support)).astype(float)
        for _ in range(4)
    )
    for attempt, repair_objective in enumerate(repair_objectives):
        repaired = linprog(
            repair_objective,
            A_eq=repair_matrix,
            b_eq=repair_rhs,
            bounds=repair_bounds,
            method="highs-ds",
        )
        if not repaired.success:
            continue
        for max_denominator in denominators:
            inequality_weights = {}
            equality_weights = {}
            for support_column, (row_type, index) in enumerate(support):
                weight = F(float(repaired.x[support_column])).limit_denominator(
                    max_denominator
                )
                if row_type == "inequality":
                    inequality_weights[index] = weight
                else:
                    equality_weights[index] = weight
            certificate = validate(
                inequality_weights,
                equality_weights,
                max_denominator,
                f"optimal_face_repair_{attempt}",
            )
            if certificate is not None:
                return certificate
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--mixed-branches", type=Path)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    mixed_branches = (
        json.loads(args.mixed_branches.read_text())[
            "selected_representations"
        ]
        if args.mixed_branches is not None
        else {}
    )
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    node_index = {node: index for index, node in enumerate(node_ids)}
    variable_count = len(node_ids) * coalition_count

    def column(node: int, coalition: int) -> int:
        return node_index[node] * coalition_count + coalition

    inequalities: list[dict[int, float]] = []
    inequality_rhs: list[float] = []
    inequality_types: list[str] = []
    inequality_metadata: list[tuple[Any, ...]] = []
    equalities: list[dict[int, float]] = []
    equality_rhs: list[float] = []
    equality_types: list[str] = []
    equality_metadata: list[tuple[Any, ...]] = []

    facets, _ = load_facets()
    for node in node_ids:
        equalities.append({column(node, 0): 1.0})
        equality_rhs.append(0.0)
        equality_types.append("empty")
        equality_metadata.append(("empty", node))
        equalities.append({column(node, grand): 1.0})
        equality_rhs.append(1.0)
        equality_types.append("grand")
        equality_metadata.append(("grand", node))
        for facet_index, facet in enumerate(facets):
            inequalities.append(
                {
                    column(node, coalition): -float(coefficient)
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                }
            )
            inequality_rhs.append(0.0)
            inequality_types.append("exact_cone")
            inequality_metadata.append(("exact_cone", node, facet_index))
        for coalition in range(coalition_count):
            for player in range(n):
                if coalition >> player & 1:
                    continue
                successor = coalition | (1 << player)
                inequalities.append(
                    {
                        column(node, coalition): 1.0,
                        column(node, successor): -1.0,
                    }
                )
                inequality_rhs.append(0.0)
                inequality_types.append("game_monotonicity")
                inequality_metadata.append(
                    ("game_monotonicity", node, coalition, successor)
                )

    protected_players: dict[tuple[int, int], set[int]] = defaultdict(set)
    for path in report["terminal_path_decomposition"]:
        protected_players[
            int(path["source"]), int(path["sink"])
        ].add(int(path["player"]))
    for (source, sink), players in protected_players.items():
        player_mask = sum(1 << player for player in players)
        for coalition in range(coalition_count):
            row = {
                column(source, coalition): 1.0,
                column(sink, coalition): -1.0,
            }
            if coalition & player_mask != player_mask:
                equalities.append(row)
                equality_rhs.append(0.0)
                equality_types.append("protected_invariance")
                equality_metadata.append(
                    ("protected_invariance", source, sink, coalition)
                )
            else:
                inequalities.append(row)
                inequality_rhs.append(0.0)
                inequality_types.append("directed_monotonicity")
                inequality_metadata.append(
                    ("directed_monotonicity", source, sink, coalition)
                )

    objective = np.zeros(variable_count)
    objective_exact = [F(0)] * variable_count
    objective_constant = F(0)
    for terminal in terminals:
        node = int(terminal["node"])
        raw_coefficient = terminal.get("coefficient")
        if raw_coefficient is not None:
            coefficient = F(raw_coefficient)
            coalition = int(terminal["coalition"])
            if coefficient > 0:
                objective[column(node, coalition)] -= float(coefficient)
                objective_exact[column(node, coalition)] -= coefficient
            elif coefficient < 0:
                mass = -coefficient
                objective[column(node, grand ^ coalition)] -= float(mass)
                objective_exact[column(node, grand ^ coalition)] -= mass
                objective_constant -= mass
        else:
            divergence = [F(value) for value in terminal["scaled_divergence"]]
            levels = sorted(set(divergence))
            if len(levels) == 2:
                low, high = levels
                coalition = sum(
                    1 << player
                    for player, value in enumerate(divergence)
                    if value == high
                )
                linear_terms = [(coalition, high - low), (grand, low)]
            else:
                selected = mixed_branches.get(str(node))
                if selected is None:
                    raise RuntimeError(
                        f"terminal {node} needs a selected mixed branch"
                    )
                linear_terms = [
                    (grand, F(selected["mu"]))
                ] + [
                    (int(term["coalition"]), F(term["weight"]))
                    for term in selected["coalitions"]
                ]
            for coalition, weight in linear_terms:
                objective[column(node, coalition)] -= float(weight)
                objective_exact[column(node, coalition)] -= weight

    result = linprog(
        objective,
        A_ub=sparse_matrix(inequalities, variable_count),
        b_ub=np.asarray(inequality_rhs),
        A_eq=sparse_matrix(equalities, variable_count),
        b_eq=np.asarray(equality_rhs),
        bounds=[(None, None)] * variable_count,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    exact_dual = reconstruct_exact_dual(
        result,
        inequalities,
        inequality_rhs,
        inequality_metadata,
        equalities,
        equality_rhs,
        equality_metadata,
        objective_exact,
        objective_constant,
    )
    maximum = -result.fun + float(objective_constant)
    nonzero_inequality_marginals = [
        float(value)
        for value in result.ineqlin.marginals
        if abs(float(value)) > 1e-9
    ]
    nonzero_equality_marginals = [
        float(value)
        for value in result.eqlin.marginals
        if abs(float(value)) > 1e-9
    ]
    active_inequality_types: dict[str, int] = defaultdict(int)
    for row_type, residual in zip(
        inequality_types, result.ineqlin.residual, strict=True
    ):
        if abs(float(residual)) <= 1e-8:
            active_inequality_types[row_type] += 1

    payload: dict[str, Any] = {
        "status": (
            "terminal_topology_relaxation_exact"
            if exact_dual is not None
            else "terminal_topology_relaxation_float"
        ),
        "source": str(args.terminal_report),
        "terminal_count": len(terminals),
        "protected_terminal_pair_count": len(protected_players),
        "variable_count": variable_count,
        "inequality_count": len(inequalities),
        "equality_count": len(equalities),
        "maximum_lower_expectation_objective_float": maximum,
        "safe_for_every_terminal_realization": (
            maximum <= 1e-9 and exact_dual is not None
        ),
        "exact_dual_verified": exact_dual is not None,
        "nonzero_inequality_dual_count": len(
            nonzero_inequality_marginals
        ),
        "nonzero_equality_dual_count": len(nonzero_equality_marginals),
        "dual_reconstruction_max_denominator_at_1e6": max(
            (
                F(value).limit_denominator(1_000_000).denominator
                for value in (
                    nonzero_inequality_marginals
                    + nonzero_equality_marginals
                )
            ),
            default=1,
        ),
        "active_inequality_types": dict(active_inequality_types),
        "mixed_terminal_representations": mixed_branches,
        "exact_dual_certificate": exact_dual,
        "interpretation": (
            "Every exact monotone terminal-game realization of the archived "
            "protected source-sink paths is feasible in this relaxation. A "
            "nonpositive maximum is therefore a topology-level no-go proof; "
            "a positive maximum means internal path constraints are essential."
        ),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
