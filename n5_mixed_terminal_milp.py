#!/usr/bin/env python3
"""Globally optimize mixed terminal lower expectations with an exact-basis MILP."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from n5_mixed_terminal_search import GameCone, Topology, four_cycle_topology


F = Fraction


def rank(matrix: list[list[F]]) -> int:
    if not matrix:
        return 0
    work = [list(row) for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                value - multiple * pivot_value
                for value, pivot_value in zip(
                    work[row], work[pivot_row], strict=True
                )
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def solve_columns(
    columns: tuple[tuple[F, ...], ...], rhs: tuple[F, ...]
) -> tuple[F, ...] | None:
    variable_count = len(columns)
    augmented = [
        [columns[column][row] for column in range(variable_count)] + [rhs[row]]
        for row in range(len(rhs))
    ]
    pivot_row = 0
    pivots: list[int] = []
    for column in range(variable_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, len(augmented))
                if augmented[row][column]
            ),
            None,
        )
        if pivot is None:
            return None
        augmented[pivot_row], augmented[pivot] = (
            augmented[pivot], augmented[pivot_row]
        )
        scale = augmented[pivot_row][column]
        augmented[pivot_row] = [
            value / scale for value in augmented[pivot_row]
        ]
        for row in range(len(augmented)):
            if row == pivot_row or not augmented[row][column]:
                continue
            multiple = augmented[row][column]
            augmented[row] = [
                value - multiple * pivot_value
                for value, pivot_value in zip(
                    augmented[row], augmented[pivot_row], strict=True
                )
            ]
        pivots.append(column)
        pivot_row += 1
    for row in range(pivot_row, len(augmented)):
        if all(not augmented[row][column] for column in range(variable_count)):
            if augmented[row][-1]:
                return None
    return tuple(augmented[row][-1] for row in range(variable_count))


def lower_dual_vertices(
    divergence: tuple[int, ...], n: int
) -> list[dict[str, Any]]:
    grand = (1 << n) - 1
    levels = sorted(set(divergence))
    if len(levels) <= 2:
        low = F(levels[0])
        high = F(levels[-1])
        coalition = sum(
            1 << player
            for player, value in enumerate(divergence)
            if value == high
        )
        coefficients = (
            {} if high == low else {coalition: high - low}
        )
        return [{"constant": low, "coefficients": coefficients}]
    coalitions = list(range(1, grand))
    columns = {
        coalition: tuple(
            F(int(coalition >> player & 1) - int(coalition & 1))
            for player in range(1, n)
        )
        for coalition in coalitions
    }
    rhs = tuple(F(divergence[player] - divergence[0]) for player in range(1, n))
    vertices: dict[tuple[tuple[int, F], ...], dict[str, Any]] = {}
    if not any(rhs):
        vertices[tuple()] = {
            "constant": F(divergence[0]),
            "coefficients": {},
        }
    for support_size in range(1, n):
        for support in itertools.combinations(coalitions, support_size):
            selected = tuple(columns[coalition] for coalition in support)
            matrix = [
                [selected[column][row] for column in range(support_size)]
                for row in range(n - 1)
            ]
            if rank(matrix) != support_size:
                continue
            solution = solve_columns(selected, rhs)
            if solution is None or any(value < 0 for value in solution):
                continue
            weights = tuple(
                (coalition, value)
                for coalition, value in zip(support, solution, strict=True)
                if value
            )
            if len(weights) != support_size:
                continue
            constant = F(divergence[0]) - sum(
                value for coalition, value in weights if coalition & 1
            )
            vertices[weights] = {
                "constant": constant,
                "coefficients": dict(weights),
            }
    return list(vertices.values())


def solve(topology: Topology, n: int) -> dict[str, Any]:
    cone = GameCone(n, topology)
    candidate_sets = [
        lower_dual_vertices(divergence, n)
        for divergence in topology.divergence
    ]
    game_variables = cone.variable_count
    z_offset = game_variables
    binary_offsets: list[int] = []
    variable_count = game_variables + topology.node_count
    for candidates in candidate_sets:
        binary_offsets.append(variable_count)
        variable_count += len(candidates)

    rows: list[dict[int, float]] = []
    lower_bounds: list[float] = []
    upper_bounds: list[float] = []
    for row, bound in zip(
        cone.inequalities, cone.inequality_rhs, strict=True
    ):
        rows.append(dict(row))
        lower_bounds.append(-np.inf)
        upper_bounds.append(bound)
    for row, bound in zip(cone.equalities, cone.equality_rhs, strict=True):
        rows.append(dict(row))
        lower_bounds.append(bound)
        upper_bounds.append(bound)

    for node, (divergence, candidates) in enumerate(
        zip(topology.divergence, candidate_sets, strict=True)
    ):
        z = z_offset + node
        minimum = float(min(divergence))
        maximum = float(max(divergence))
        for candidate_index, candidate in enumerate(candidates):
            constant = float(candidate["constant"])
            coefficients = candidate["coefficients"]
            # affine(v) - z <= 0, hence z is at least every dual value.
            row = {z: -1.0}
            for coalition, coefficient in coefficients.items():
                row[cone.column(node, coalition)] = float(coefficient)
            rows.append(row)
            lower_bounds.append(-np.inf)
            upper_bounds.append(-constant)

            affine_minimum = constant
            big_m = max(1.0, maximum - affine_minimum)
            # z - affine(v) + M*b <= M.  When b=1 this selects the
            # maximizing affine piece; otherwise the row is relaxed.
            row = {
                z: 1.0,
                binary_offsets[node] + candidate_index: big_m,
            }
            for coalition, coefficient in coefficients.items():
                row[cone.column(node, coalition)] = -float(coefficient)
            rows.append(row)
            lower_bounds.append(-np.inf)
            upper_bounds.append(big_m + constant)
        rows.append(
            {
                binary_offsets[node] + candidate_index: 1.0
                for candidate_index in range(len(candidates))
            }
        )
        lower_bounds.append(1.0)
        upper_bounds.append(1.0)

    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[float] = []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                row_indices.append(row_index)
                column_indices.append(column)
                values.append(value)
    matrix = coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(rows), variable_count),
    ).tocsr()
    objective = np.zeros(variable_count)
    objective[z_offset : z_offset + topology.node_count] = -1.0
    variable_lower = np.full(variable_count, -np.inf)
    variable_upper = np.full(variable_count, np.inf)
    variable_lower[:game_variables] = 0.0
    variable_upper[:game_variables] = 1.0
    for node, divergence in enumerate(topology.divergence):
        variable_lower[z_offset + node] = min(divergence)
        variable_upper[z_offset + node] = max(divergence)
        start = binary_offsets[node]
        stop = start + len(candidate_sets[node])
        variable_lower[start:stop] = 0.0
        variable_upper[start:stop] = 1.0
    integrality = np.zeros(variable_count)
    for node, candidates in enumerate(candidate_sets):
        start = binary_offsets[node]
        integrality[start : start + len(candidates)] = 1
    result = milp(
        objective,
        integrality=integrality,
        bounds=Bounds(variable_lower, variable_upper),
        constraints=LinearConstraint(
            matrix, np.asarray(lower_bounds), np.asarray(upper_bounds)
        ),
        options={"time_limit": 600.0, "mip_rel_gap": 0.0},
    )
    payload: dict[str, Any] = {
        "status": "optimal" if result.success else "incomplete",
        "message": result.message,
        "objective": None if result.fun is None else -float(result.fun),
        "candidate_counts": [len(candidates) for candidates in candidate_sets],
        "node_count": topology.node_count,
        "arcs": [list(arc) for arc in topology.arcs],
        "divergence": [list(vector) for vector in topology.divergence],
        "mip_gap": getattr(result, "mip_gap", None),
        "mip_node_count": getattr(result, "mip_node_count", None),
    }
    if result.x is not None:
        payload["node_lower_expectations"] = [
            float(result.x[z_offset + node])
            for node in range(topology.node_count)
        ]
        payload["games"] = [
            [
                float(result.x[cone.column(node, coalition)])
                for coalition in range(cone.coalition_count)
            ]
            for node in range(topology.node_count)
        ]
    return payload


def complementary_mixed_fork_topology(n: int) -> Topology:
    if n < 3:
        raise ValueError("the complementary mixed fork needs at least three players")
    first = (1, 2) + (0,) * (n - 2)
    second = (1, 0) + (2,) * (n - 2)
    root = tuple(left + right for left, right in zip(first, second, strict=True))
    first_support = sum(1 << player for player, value in enumerate(first) if value)
    second_support = sum(1 << player for player, value in enumerate(second) if value)
    return Topology(
        3,
        ((0, 1, first_support), (0, 2, second_support)),
        (root, tuple(-value for value in first), tuple(-value for value in second)),
    )


def indicator_hourglass_topology(masks: tuple[int, int, int, int], n: int) -> Topology:
    grand = (1 << n) - 1
    if any(mask <= 0 or mask >= grand for mask in masks):
        raise ValueError("hourglass masks must be nonempty proper coalitions")
    incoming_left, incoming_right, outgoing_left, outgoing_right = masks
    divergence = [[0] * n for _ in range(5)]
    arcs = (
        (0, 2, incoming_left),
        (1, 2, incoming_right),
        (2, 3, outgoing_left),
        (2, 4, outgoing_right),
    )
    for player in range(n):
        if incoming_left >> player & 1:
            divergence[0][player] += 1
            divergence[2][player] -= 1
        if incoming_right >> player & 1:
            divergence[1][player] += 1
            divergence[2][player] -= 1
        if outgoing_left >> player & 1:
            divergence[2][player] += 1
            divergence[3][player] -= 1
        if outgoing_right >> player & 1:
            divergence[2][player] += 1
            divergence[4][player] -= 1
    return Topology(5, arcs, tuple(tuple(vector) for vector in divergence))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", type=int, nargs=4, default=(1, 2, 2, 1))
    parser.add_argument("--complementary-mixed-fork", action="store_true")
    parser.add_argument("--hourglass-masks", type=int, nargs=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    topology = (
        indicator_hourglass_topology(tuple(args.hourglass_masks), 5)
        if args.hourglass_masks is not None
        else complementary_mixed_fork_topology(5)
        if args.complementary_mixed_fork
        else four_cycle_topology(tuple(args.weights), 5)
    )
    payload = solve(topology, 5)
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0 if payload["status"] == "optimal" else 1


if __name__ == "__main__":
    raise SystemExit(main())
