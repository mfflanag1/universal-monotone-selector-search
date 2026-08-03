#!/usr/bin/env python3
"""Prove a fixed-grand K2,3 terminal family has a common core.

The common-core envelope is balanced iff every extreme balanced collection
has worth at most the common grand worth.  Only sink coordinates not forced
equal by the two protected source rows require a choice of which sink attains
the envelope.  The script exhausts those finite choices and rationally audits
the LP dual for every case.
"""

from __future__ import annotations

import argparse
import itertools
import json
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize._highspy._core import (
    HighsLp,
    MatrixFormat,
    ObjSense,
    _Highs,
    kHighsInf,
)
from scipy.sparse import vstack

from n5_mixed_terminal_search import (
    GameCone,
    complete_bipartite_rectangular_topology,
)


F = Fraction
N = 5
GRAND = (1 << N) - 1


DETERMINANT_TERMS = {
    size: [
        (
            permutation,
            -1
            if sum(
                permutation[left] > permutation[right]
                for left in range(size)
                for right in range(left + 1, size)
            )
            % 2
            else 1,
        )
        for permutation in itertools.permutations(range(size))
    ]
    for size in range(1, N + 1)
}


def integer_determinant(matrix: list[list[int]]) -> int:
    total = 0
    for permutation, sign in DETERMINANT_TERMS[len(matrix)]:
        product = 1
        for row, column in enumerate(permutation):
            product *= matrix[row][column]
            if not product:
                break
        total += sign * product
    return total


def balanced_vertices() -> list[tuple[tuple[int, F], ...]]:
    columns = {
        coalition: tuple(
            int(coalition >> player & 1) for player in range(N)
        )
        for coalition in range(1, GRAND)
    }
    vertices: set[tuple[tuple[int, F], ...]] = set()
    for support_size in range(1, N + 1):
        for support in itertools.combinations(range(1, GRAND), support_size):
            exact: tuple[F, ...] | None = None
            for selected_rows in itertools.combinations(range(N), support_size):
                square = [
                    [columns[coalition][row] for coalition in support]
                    for row in selected_rows
                ]
                denominator = integer_determinant(square)
                if not denominator:
                    continue
                numerators = []
                for replaced_column in range(support_size):
                    replaced = [list(row) for row in square]
                    for row in replaced:
                        row[replaced_column] = 1
                    numerators.append(integer_determinant(replaced))
                exact = tuple(
                    F(numerator, denominator) for numerator in numerators
                )
                break
            if exact is None or any(value <= 0 for value in exact):
                continue
            if not all(
                sum(
                    exact[index] * int(coalition >> player & 1)
                    for index, coalition in enumerate(support)
                )
                == 1
                for player in range(N)
            ):
                continue
            vertices.add(tuple(zip(support, exact, strict=True)))
    return sorted(vertices, key=str)


def load_or_build_balanced_vertices(
    cache: Path | None,
) -> list[tuple[tuple[int, F], ...]]:
    if cache is not None and cache.exists():
        payload = json.loads(cache.read_text())
        vertices = [
            tuple((int(coalition), F(weight)) for coalition, weight in vertex)
            for vertex in payload["vertices"]
        ]
        if payload.get("status") != "n5_balanced_vertices_exact":
            raise RuntimeError("balanced-vertex cache has the wrong status")
        if len(vertices) != 1_291:
            raise RuntimeError("balanced-vertex cache has the wrong count")
        for vertex in vertices:
            if any(weight <= 0 for _, weight in vertex):
                raise RuntimeError("balanced-vertex cache has a nonpositive weight")
            if any(
                sum(
                    weight * int(coalition >> player & 1)
                    for coalition, weight in vertex
                )
                != 1
                for player in range(N)
            ):
                raise RuntimeError("balanced-vertex cache has an invalid vertex")
        return vertices
    vertices = balanced_vertices()
    if cache is not None:
        cache.write_text(
            json.dumps(
                {
                    "status": "n5_balanced_vertices_exact",
                    "enumeration_method": "integer_determinants_and_cramers_rule",
                    "vertex_count": len(vertices),
                    "vertices": [
                        [[coalition, str(weight)] for coalition, weight in vertex]
                        for vertex in vertices
                    ],
                },
                indent=2,
            )
            + "\n"
        )
    return vertices


def varying_sink_coordinates(masks: tuple[int, ...]) -> list[int]:
    rows = masks[:3], masks[3:]
    varying = []
    for coalition in range(1, GRAND):
        equality_graph = [set() for _ in range(3)]
        for left in range(3):
            for right in range(left + 1, 3):
                if any(
                    coalition & row[left] != row[left]
                    and coalition & row[right] != row[right]
                    for row in rows
                ):
                    equality_graph[left].add(right)
                    equality_graph[right].add(left)
        seen = {0}
        stack = [0]
        while stack:
            node = stack.pop()
            for neighbor in equality_graph[node]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        if len(seen) != 3:
            varying.append(coalition)
    return varying


def exact_dual_denominator(
    row_dual: np.ndarray,
    exact_rows: list[dict[int, F]],
    exact_rhs: list[F],
    inequality_count: int,
    exact_cost: list[F],
) -> int | None:
    support = [
        index for index, value in enumerate(row_dual) if abs(float(value)) > 1e-9
    ]
    for denominator in (
        100,
        1_000,
        10_000,
        100_000,
        1_000_000,
        100_000_000,
    ):
        weights = {
            index: F(float(row_dual[index])).limit_denominator(denominator)
            for index in support
        }
        if any(
            weight > 0
            for index, weight in weights.items()
            if index < inequality_count
        ):
            continue
        stationarity = [F(0)] * len(exact_cost)
        dual_objective = F(0)
        for index, weight in weights.items():
            for column, coefficient in exact_rows[index].items():
                stationarity[column] += weight * coefficient
            dual_objective += weight * exact_rhs[index]
        if stationarity == exact_cost and dual_objective == -1:
            return denominator
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--masks", type=int, nargs=6, default=(4, 21, 7, 13, 1, 11)
    )
    parser.add_argument("--balanced-cache", type=Path)
    parser.add_argument("--balanced-cache-only", action="store_true")
    parser.add_argument("--progress-every", type=int, default=500)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    masks = tuple(args.masks)
    topology = complete_bipartite_rectangular_topology(
        masks, (1,) * 6, 2, 3, N
    )
    cone = GameCone(N, topology)
    vertices = load_or_build_balanced_vertices(args.balanced_cache)
    if args.balanced_cache_only:
        print(
            json.dumps(
                {
                    "status": "n5_balanced_vertices_exact",
                    "vertex_count": len(vertices),
                    "cache": (
                        None
                        if args.balanced_cache is None
                        else str(args.balanced_cache)
                    ),
                },
                indent=2,
            )
        )
        return 0
    varying = set(varying_sink_coordinates(masks))

    matrix = vstack((cone.a_ub, cone.a_eq)).tocsc()
    exact_rows = [
        {column: F(value).limit_denominator() for column, value in row.items()}
        for row in cone.inequalities + cone.equalities
    ]
    exact_rhs = [
        F(value).limit_denominator()
        for value in cone.inequality_rhs + cone.equality_rhs
    ]
    lp = HighsLp()
    lp.num_col_ = cone.variable_count
    lp.num_row_ = matrix.shape[0]
    lp.col_cost_ = np.zeros(cone.variable_count)
    lp.col_lower_ = np.full(cone.variable_count, -kHighsInf)
    lp.col_upper_ = np.full(cone.variable_count, kHighsInf)
    lp.row_lower_ = np.concatenate(
        (np.full(cone.a_ub.shape[0], -kHighsInf), cone.b_eq)
    )
    lp.row_upper_ = np.concatenate((cone.b_ub, cone.b_eq))
    lp.a_matrix_.format_ = MatrixFormat.kColwise
    lp.a_matrix_.num_col_ = cone.variable_count
    lp.a_matrix_.num_row_ = matrix.shape[0]
    lp.a_matrix_.start_ = matrix.indptr.astype(np.int32)
    lp.a_matrix_.index_ = matrix.indices.astype(np.int32)
    lp.a_matrix_.value_ = matrix.data
    lp.sense_ = ObjSense.kMinimize
    highs = _Highs()
    highs.setOptionValue("output_flag", False)
    highs.setOptionValue("solver", "simplex")
    highs.setOptionValue("simplex_strategy", 1)
    if highs.passModel(lp).name != "kOk" or highs.run().name != "kOk":
        raise RuntimeError("HiGHS failed to initialize the common-core sweep")

    indices = np.arange(cone.variable_count, dtype=np.int32)
    total = sum(
        3 ** sum(coalition in varying for coalition, _ in vertex)
        for vertex in vertices
    )
    solved = 0
    maximum = -float("inf")
    failures: list[dict[str, object]] = []
    maximum_denominator = 1
    started = time.monotonic()
    for vertex_index, vertex in enumerate(vertices):
        variable_terms = [coalition for coalition, _ in vertex if coalition in varying]
        for selected_sinks in itertools.product(range(3), repeat=len(variable_terms)):
            selection = dict(zip(variable_terms, selected_sinks, strict=True))
            exact_cost = [F(0)] * cone.variable_count
            for coalition, weight in vertex:
                sink = selection.get(coalition, 0)
                exact_cost[cone.column(2 + sink, coalition)] -= weight
            cost = np.asarray([float(value) for value in exact_cost])
            if highs.changeColsCost(cone.variable_count, indices, cost).name != "kOk":
                raise RuntimeError("HiGHS rejected an objective update")
            if highs.run().name != "kOk":
                raise RuntimeError(
                    f"HiGHS failed at balanced vertex {vertex_index}"
                )
            value = -float(highs.getObjectiveValue())
            maximum = max(maximum, value)
            denominator = exact_dual_denominator(
                np.asarray(highs.getSolution().row_dual),
                exact_rows,
                exact_rhs,
                cone.a_ub.shape[0],
                exact_cost,
            )
            if denominator is None:
                failures.append(
                    {
                        "balanced_vertex": vertex_index,
                        "selected_sinks": list(selected_sinks),
                        "maximum_float": value,
                    }
                )
            else:
                maximum_denominator = max(maximum_denominator, denominator)
            solved += 1
            if args.progress_every and solved % args.progress_every == 0:
                print(
                    json.dumps(
                        {
                            "solved": solved,
                            "total": total,
                            "maximum": maximum,
                            "exact_failures": len(failures),
                        }
                    ),
                    flush=True,
                )
    payload = {
        "status": (
            "k23_common_core_exact"
            if solved == total and not failures
            else "k23_common_core_float"
        ),
        "masks": list(masks),
        "balanced_vertex_count": len(vertices),
        "balanced_vertex_enumeration": "integer_determinants_and_cramers_rule",
        "varying_sink_coordinates": sorted(varying),
        "envelope_assignment_count": total,
        "solved": solved,
        "maximum_balanced_worth_float": maximum,
        "exact_certificate_count": solved - len(failures),
        "exact_failure_count": len(failures),
        "first_exact_failures": failures[:20],
        "exact_reconstruction_max_denominator": maximum_denominator,
        "common_core_gap_exact": "0" if solved == total and not failures else None,
        "elapsed_seconds": time.monotonic() - started,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0 if payload["status"] == "k23_common_core_exact" else 1


if __name__ == "__main__":
    raise SystemExit(main())
