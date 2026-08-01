#!/usr/bin/env python3
"""Exactify an envelope-constrained n=5 topology from its reduced game cone."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.linalg import qr
from scipy.optimize import linprog
from sympy import Rational
from sympy.polys.matrices import DomainMatrix

from n5_facet_search import exact_archive, load_facets, margin_dual


F = Fraction
N = 5
GRAND = 31


def parse_ints(text: str) -> tuple[int, ...]:
    return tuple(int(value) for value in text.split(","))


def topology_offsets(
    raw: dict[str, Any],
) -> tuple[list[list[F]], list[tuple[int, int, int]]]:
    game_count, edges, deltas = archive_edges(raw)
    return propagate_offsets(game_count, edges, deltas)


def archive_edges(
    raw: dict[str, Any],
) -> tuple[
    int, list[tuple[int, int, int]], list[F]
]:
    normalizer = F(raw["games"][0][GRAND])
    edges = [
        (
            int(edge["lower"]),
            int(edge["upper"]),
            int(edge["coalition"]),
        )
        for edge in raw["edges"]
    ]
    deltas = [
        F(edge["delta"]) / normalizer for edge in raw["edges"]
    ]
    return len(raw["games"]), edges, deltas


def glued_archive_offsets(
    base_raw: dict[str, Any],
    gadget_raw: dict[str, Any],
    attachment: int,
    gadget_anchor: int,
    gadget_scale: F,
    gadget_permutation: tuple[int, ...],
) -> tuple[list[list[F]], list[tuple[int, int, int]]]:
    base_count, base_edges, base_deltas = archive_edges(base_raw)
    gadget_count, gadget_edges, gadget_deltas = archive_edges(
        gadget_raw
    )
    remap = {gadget_anchor: attachment}
    next_node = base_count
    for node in range(gadget_count):
        if node != gadget_anchor:
            remap[node] = next_node
            next_node += 1
    edges = list(base_edges)
    deltas = list(base_deltas)
    for (lower, upper, coalition), delta in zip(
        gadget_edges, gadget_deltas, strict=True
    ):
        edges.append(
            (
                remap[lower],
                remap[upper],
                permute_mask(coalition, gadget_permutation),
            )
        )
        deltas.append(delta * gadget_scale)
    return propagate_offsets(next_node, edges, deltas)


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(
        1 << permutation[player]
        for player in range(len(permutation))
        if mask >> player & 1
    )


def record_offsets(
    record: dict[str, Any], max_denominator: int
) -> tuple[list[list[F]], list[tuple[int, int, int]]]:
    edges = [
        tuple(int(value) for value in edge)
        for edge in record["edges"]
    ]
    game_count = 1 + max(
        max(lower, upper) for lower, upper, _ in edges
    )
    deltas = [
        F(value).limit_denominator(max_denominator)
        for value in record["bumps"]
    ]
    return propagate_offsets(game_count, edges, deltas)


def propagate_offsets(
    game_count: int,
    edges: Sequence[tuple[int, int, int]],
    deltas: Sequence[F],
) -> tuple[list[list[F]], list[tuple[int, int, int]]]:
    offsets: list[list[F] | None] = [None] * game_count
    offsets[0] = [F(0)] * (GRAND + 1)
    for _ in range(game_count):
        changed = False
        for (lower, upper, coalition), delta in zip(
            edges, deltas, strict=True
        ):
            if offsets[lower] is not None and offsets[upper] is None:
                offsets[upper] = list(offsets[lower])
                offsets[upper][coalition] += delta
                changed = True
            elif offsets[upper] is not None and offsets[lower] is None:
                offsets[lower] = list(offsets[upper])
                offsets[lower][coalition] -= delta
                changed = True
            elif offsets[lower] is not None and offsets[upper] is not None:
                expected = list(offsets[lower])
                expected[coalition] += delta
                if expected != offsets[upper]:
                    raise ValueError("inconsistent exact bump cycle")
        if not changed:
            break
    if any(offset is None for offset in offsets):
        raise ValueError("topology is not connected")
    return [offset for offset in offsets if offset is not None], edges


def reduced_constraints(
    offsets: Sequence[Sequence[F]],
    envelope_nodes: Sequence[int],
    envelope_coalitions: Sequence[int],
    envelope_gap: F,
    envelope_multiplier: F,
) -> tuple[list[list[F]], list[F]]:
    facets, _ = load_facets()
    constraints: dict[tuple[F, ...], F] = {}

    def add(row: list[F], bound: F) -> None:
        key = tuple(row)
        constraints[key] = min(bound, constraints.get(key, bound))

    for offset in offsets:
        for facet in facets:
            row = [F(0)] * GRAND
            bound = F(0)
            for coalition, coefficient in enumerate(facet):
                if coalition and coefficient:
                    row[coalition - 1] -= coefficient
                    bound += coefficient * offset[coalition]
            add(row, bound)
        for coalition in range(GRAND):
            for player in range(N):
                if coalition >> player & 1:
                    continue
                upper = coalition | (1 << player)
                row = [F(0)] * GRAND
                if coalition:
                    row[coalition - 1] += 1
                row[upper - 1] -= 1
                bound = offset[upper] - (
                    offset[coalition] if coalition else F(0)
                )
                add(row, bound)

    envelope_row = [F(0)] * GRAND
    envelope_rhs = -envelope_multiplier * (1 + envelope_gap)
    for coalition, node in zip(
        envelope_coalitions, envelope_nodes, strict=True
    ):
        envelope_row[coalition - 1] -= 1
        envelope_rhs += offsets[node][coalition]
    add(envelope_row, envelope_rhs)
    return (
        [list(row) for row in constraints],
        list(constraints.values()),
    )


def flow_objective(
    games: Sequence[Sequence[float]],
    edges: Sequence[tuple[int, int, int]],
    max_denominator: int,
) -> list[F]:
    _, inequality_weights, equality_weights, metadata = margin_dual(
        games, edges
    )
    objective = [F(0)] * GRAND
    for weight, tag in zip(
        inequality_weights, metadata, strict=True
    ):
        if weight and tag[0] == "core":
            objective[int(tag[2]) - 1] += F(
                weight
            ).limit_denominator(max_denominator)
    for weight in equality_weights:
        objective[GRAND - 1] -= F(weight).limit_denominator(
            max_denominator
        )
    return objective


def domain_matrix(
    rows: Sequence[Sequence[F]],
) -> DomainMatrix:
    values: dict[int, dict[int, Any]] = {}
    for row_index, row in enumerate(rows):
        for column_index, value in enumerate(row):
            if value:
                values.setdefault(row_index, {})[
                    column_index
                ] = Rational(value.numerator, value.denominator)
    return DomainMatrix.from_dict_sympy(
        len(rows), len(rows[0]), values
    ).to_field()


def solve_basis(
    basis_indices: Sequence[int],
    rows: Sequence[Sequence[F]],
    rhs: Sequence[F],
    objective: Sequence[F],
) -> tuple[list[F], list[F], list[F], list[F]] | None:
    equality = [F(0)] * (GRAND - 1) + [F(1)]
    basis_rows = [equality] + [
        list(rows[index]) for index in basis_indices
    ]
    basis_rhs = [F(1)] + [rhs[index] for index in basis_indices]
    if (
        np.linalg.matrix_rank(
            np.asarray(
                [
                    [float(value) for value in row]
                    for row in basis_rows
                ]
            )
        )
        < GRAND
    ):
        return None
    basis = domain_matrix(basis_rows)
    try:
        numerator, denominator = basis.solve_den(
            domain_matrix([[value] for value in basis_rhs]),
            method="rref",
        )
        dual_numerator, dual_denominator = (
            basis.transpose().solve_den(
                domain_matrix([[value] for value in objective]),
                method="rref",
            )
        )
    except Exception:
        return None
    primal_matrix = numerator.to_Matrix()
    primal_scale = F(str(denominator))
    point = [
        F(str(primal_matrix[index, 0])) / primal_scale
        for index in range(GRAND)
    ]
    dual_matrix = dual_numerator.to_Matrix()
    dual_scale = F(str(dual_denominator))
    dual = [
        F(str(dual_matrix[index, 0])) / dual_scale
        for index in range(GRAND)
    ]
    if any(weight > 0 for weight in dual[1:]):
        return None
    violations = [
        sum(
            (
                coefficient * value
                for coefficient, value in zip(
                    row, point, strict=True
                )
            ),
            start=F(0),
        )
        - bound
        for row, bound in zip(rows, rhs, strict=True)
    ]
    return point, dual, basis_rhs, violations


def exact_vertex(
    rows: Sequence[Sequence[F]],
    rhs: Sequence[F],
    objective: Sequence[F],
) -> tuple[list[F], int]:
    equality = [F(0)] * (GRAND - 1) + [F(1)]
    float_rows = np.asarray(
        [[float(value) for value in row] for row in rows]
    )
    result = linprog(
        [float(value) for value in objective],
        A_ub=float_rows,
        b_ub=[float(value) for value in rhs],
        A_eq=[[float(value) for value in equality]],
        b_eq=[1.0],
        bounds=[(None, None)] * GRAND,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    active = np.flatnonzero(result.ineqlin.residual < 1e-8)
    support = np.flatnonzero(
        abs(result.ineqlin.marginals) > 1e-10
    )
    support_set = {int(index) for index in support}
    base = np.vstack(
        [
            np.asarray(
                [[float(value) for value in equality]]
            ),
            float_rows[support],
        ]
    )
    base_rank = int(np.linalg.matrix_rank(base))
    if base_rank != len(base):
        raise RuntimeError("dependent numerical dual support")
    q_matrix, _ = qr(base.T, mode="full")
    remaining = np.asarray(
        [
            int(index)
            for index in active
            if int(index) not in support_set
        ]
    )
    needed = GRAND - base_rank
    projected = float_rows[remaining] @ q_matrix[:, base_rank:]
    _, _, pivots = qr(
        projected.T, pivoting=True, mode="economic"
    )
    basis_indices = [
        *[int(index) for index in support],
        *[
            int(index)
            for index in remaining[pivots[:needed]]
        ],
    ]

    pivot_count = 0
    while True:
        solved = solve_basis(
            basis_indices, rows, rhs, objective
        )
        if solved is None:
            raise RuntimeError("failed to reconstruct exact basis")
        point, _, _, violations = solved
        maximum_violation = max(violations)
        if maximum_violation <= 0:
            return point, pivot_count
        entering = violations.index(maximum_violation)
        candidates = []
        for leaving in range(len(basis_indices)):
            trial = list(basis_indices)
            trial[leaving] = entering
            if len(set(trial)) != len(trial):
                continue
            candidate = solve_basis(
                trial, rows, rhs, objective
            )
            if candidate is not None:
                candidates.append(
                    (max(candidate[3]), trial)
                )
        if not candidates:
            raise RuntimeError("exact dual-simplex pivot failed")
        _, basis_indices = min(
            candidates, key=lambda candidate: candidate[0]
        )
        pivot_count += 1
        if pivot_count >= 100:
            raise RuntimeError("exact dual-simplex pivot limit")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--topology", type=Path)
    parser.add_argument("--base", type=Path)
    parser.add_argument("--gadget", type=Path)
    parser.add_argument("--attachment", type=int)
    parser.add_argument("--gadget-anchor", type=int)
    parser.add_argument("--gadget-scale")
    parser.add_argument(
        "--gadget-permutation", default="0,1,2,3,4"
    )
    parser.add_argument("--float-record", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--envelope-nodes", required=True)
    parser.add_argument(
        "--envelope-coalitions", default="30,29,27,23,15"
    )
    parser.add_argument("--envelope-gap", required=True)
    parser.add_argument(
        "--envelope-multiplier", default="4"
    )
    parser.add_argument(
        "--max-denominator", type=int, default=1_000_000
    )
    args = parser.parse_args()

    captured = json.loads(args.float_record.read_text())
    record = captured.get("record", captured)
    if args.base is not None or args.gadget is not None:
        if None in (
            args.base,
            args.gadget,
            args.attachment,
            args.gadget_anchor,
            args.gadget_scale,
        ):
            raise ValueError(
                "exact glue requires base, gadget, attachment, "
                "gadget anchor, and gadget scale"
            )
        base = json.loads(args.base.read_text())
        gadget = json.loads(args.gadget.read_text())
        gadget_permutation = parse_ints(
            args.gadget_permutation
        )
        if sorted(gadget_permutation) != list(range(N)):
            raise ValueError(
                "gadget permutation must contain 0,1,2,3,4"
            )
        offsets, edges = glued_archive_offsets(
            base["family"],
            gadget["family"],
            args.attachment,
            args.gadget_anchor,
            F(args.gadget_scale),
            gadget_permutation,
        )
    else:
        if args.topology is None:
            raise ValueError("topology or exact glue inputs required")
        topology = json.loads(args.topology.read_text())
        if topology.get("family") is not None:
            offsets, edges = topology_offsets(topology["family"])
        else:
            topology_record = topology.get("record", topology)
            offsets, edges = record_offsets(
                topology_record, args.max_denominator
            )
    record_edges = [
        tuple(int(value) for value in edge)
        for edge in record["edges"]
    ]
    if record_edges != edges:
        mismatch_count = sum(
            left != right
            for left, right in zip(record_edges, edges)
        ) + abs(len(record_edges) - len(edges))
        first_mismatch = next(
            (
                (index, left, right)
                for index, (left, right) in enumerate(
                    zip(record_edges, edges)
                )
                if left != right
            ),
            None,
        )
        raise ValueError(
            "float record topology differs from exact topology: "
            f"{mismatch_count} mismatches; first={first_mismatch}"
        )
    if len(record["best_games_float"]) != len(offsets):
        raise ValueError(
            "float record game count differs from exact topology"
        )
    envelope_nodes = parse_ints(args.envelope_nodes)
    envelope_coalitions = parse_ints(
        args.envelope_coalitions
    )
    if len(envelope_nodes) != len(envelope_coalitions):
        raise ValueError("envelope node and coalition counts differ")
    rows, rhs = reduced_constraints(
        offsets,
        envelope_nodes,
        envelope_coalitions,
        F(args.envelope_gap),
        F(args.envelope_multiplier),
    )
    objective = flow_objective(
        record["best_games_float"],
        edges,
        args.max_denominator,
    )
    base, pivot_count = exact_vertex(rows, rhs, objective)
    exact_games = [
        [
            str(
                (F(0) if coalition == 0 else base[coalition - 1])
                + offset[coalition]
            )
            for coalition in range(GRAND + 1)
        ]
        for offset in offsets
    ]
    record["best_games_exact"] = exact_games
    archive = exact_archive(
        record,
        "n5_exactified_envelope_constrained_topology",
        args.max_denominator,
    )
    archive["exact_game_vertex_reconstruction"] = {
        "topology": str(args.topology) if args.topology else None,
        "base": str(args.base) if args.base else None,
        "gadget": str(args.gadget) if args.gadget else None,
        "attachment": args.attachment,
        "gadget_anchor": args.gadget_anchor,
        "gadget_scale": (
            str(F(args.gadget_scale))
            if args.gadget_scale is not None
            else None
        ),
        "gadget_permutation": (
            list(gadget_permutation)
            if args.base is not None
            else None
        ),
        "float_record": str(args.float_record),
        "reduced_constraint_count": len(rows),
        "exact_dual_simplex_pivots": pivot_count,
        "envelope_nodes": list(envelope_nodes),
        "envelope_coalitions": list(envelope_coalitions),
        "envelope_gap": str(F(args.envelope_gap)),
        "envelope_multiplier": str(F(args.envelope_multiplier)),
        "maximum_base_denominator": max(
            value.denominator for value in base
        ),
    }
    args.output.write_text(json.dumps(archive, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: archive[key]
                for key in (
                    "status",
                    "node_count",
                    "edge_count",
                    "common_core_budget_gap_exact",
                    "max_min_monotonicity_margin_exact",
                    "box_max_min_monotonicity_margin_exact",
                    "non_atomic_facet_tax_exact",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
