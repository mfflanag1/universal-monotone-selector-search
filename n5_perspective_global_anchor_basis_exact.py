#!/usr/bin/env python3
"""Solve a saved global perspective support by sparse exact elimination."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import qr
from scipy.sparse import coo_matrix, csr_matrix, vstack
from sympy import Matrix, Rational, SparseMatrix, linsolve

from n5_facet_search import load_facets


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("anchor_search", type=Path)
    parser.add_argument("float_support", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--diagnostics-only", action="store_true")
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    anchor = json.loads(args.anchor_search.read_text())[
        "best_complete_anchor"
    ]
    support = json.loads(args.float_support.read_text())
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    node_index = {node: index for index, node in enumerate(node_ids)}
    variable_count = len(node_ids) * coalition_count
    mixed_nodes = {int(node) for node in anchor["links"]}
    facets, _ = load_facets()

    def column(node: int, coalition: int) -> int:
        return node_index[node] * coalition_count + coalition

    protected_players: dict[tuple[int, int], set[int]] = defaultdict(set)
    for path in report["terminal_path_decomposition"]:
        protected_players[
            int(path["source"]), int(path["sink"])
        ].add(int(path["player"]))

    def row_from_metadata(metadata: list[Any]) -> tuple[dict[int, int], int]:
        row_type = metadata[0]
        if row_type == "exact_cone":
            _, node, facet_index = metadata
            facet = facets[int(facet_index)]
            return (
                {
                    column(int(node), coalition): -int(coefficient)
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                },
                0,
            )
        if row_type == "game_monotonicity":
            _, node, coalition, successor = metadata
            return (
                {
                    column(int(node), int(coalition)): 1,
                    column(int(node), int(successor)): -1,
                },
                0,
            )
        if row_type == "directed_monotonicity":
            _, source, sink, coalition = metadata
            return (
                {
                    column(int(source), int(coalition)): 1,
                    column(int(sink), int(coalition)): -1,
                },
                0,
            )
        if row_type == "empty":
            _, node = metadata
            return ({column(int(node), 0): 1}, 0)
        if row_type == "grand":
            _, node = metadata
            return ({column(int(node), grand): 1}, 1)
        if row_type == "protected_invariance":
            _, source, sink, coalition = metadata
            return (
                {
                    column(int(source), int(coalition)): 1,
                    column(int(sink), int(coalition)): -1,
                },
                0,
            )
        raise RuntimeError(f"unsupported support row {row_type}")

    entries = support["inequalities"] + support["equalities"]
    inequality_count = len(support["inequalities"])
    rows = [row_from_metadata(active["row"])[0] for active in entries]
    rhs = [row_from_metadata(active["row"])[1] for active in entries]
    row_indices = []
    columns = []
    values = []
    for support_column, row in enumerate(rows):
        for variable, coefficient in row.items():
            row_indices.append(variable)
            columns.append(support_column)
            values.append(float(coefficient))
    stationarity = coo_matrix(
        (values, (row_indices, columns)),
        shape=(variable_count, len(rows)),
    ).tocsr()

    objective = [F(0)] * variable_count
    for terminal in terminals:
        node = int(terminal["node"])
        if node in mixed_nodes:
            continue
        raw_coefficient = terminal.get("coefficient")
        if raw_coefficient is not None:
            coefficient = F(raw_coefficient)
            coalition = int(terminal["coalition"])
            if coefficient > 0:
                objective[column(node, coalition)] -= coefficient
            elif coefficient < 0:
                objective[column(node, grand ^ coalition)] -= -coefficient
        else:
            divergence = [F(value) for value in terminal["scaled_divergence"]]
            low, high = sorted(set(divergence))
            coalition = sum(
                1 << player
                for player, value in enumerate(divergence)
                if value == high
            )
            objective[column(node, coalition)] -= high - low
            objective[column(node, grand)] -= low
    links = {
        int(node): [F(value) for value in link]
        for node, link in anchor["links"].items()
    }
    for node, link in links.items():
        for coalition, value in enumerate(link):
            objective[column(node, coalition)] -= value
    objective_constant = sum(
        (
            F(terminal["coefficient"])
            for terminal in terminals
            if terminal.get("coefficient") is not None
            and F(terminal["coefficient"]) < 0
        ),
        F(0),
    )
    simplex = {
        int(node): F(value)
        for node, value in anchor["simplex"].items()
    }
    target_dual_objective = objective_constant - sum(
        simplex.values(), F(0)
    )
    augmented_stationarity = vstack(
        (
            stationarity,
            csr_matrix(np.asarray(rhs, dtype=float).reshape(1, -1)),
        ),
        format="csr",
    )
    augmented_objective = objective + [target_dual_objective]

    _, triangular, pivots = qr(
        augmented_stationarity.transpose().toarray(),
        mode="economic",
        pivoting=True,
        check_finite=False,
    )
    diagonal = np.abs(np.diag(triangular))
    tolerance = (
        max(augmented_stationarity.shape)
        * np.finfo(float).eps
        * (diagonal[0] if len(diagonal) else 0.0)
    )
    rank = int(np.count_nonzero(diagonal > tolerance))
    diagnostics = {
        "equation_count": len(augmented_objective),
        "support_variable_count": len(rows),
        "numeric_rank": rank,
        "rank_tolerance": tolerance,
        "smallest_selected_diagonal": (
            float(diagonal[rank - 1]) if rank else None
        ),
    }
    print(json.dumps(diagnostics, indent=2), flush=True)
    if args.diagnostics_only:
        return 0
    if rank != len(rows):
        raise RuntimeError("saved global support does not have full column rank")

    selected_equations = [int(value) for value in pivots[:rank]]
    square = augmented_stationarity[selected_equations, :].tocoo()
    exact_square = SparseMatrix(
        rank,
        rank,
        {
            (int(row), int(column)): int(round(value))
            for row, column, value in zip(
                square.row, square.col, square.data, strict=True
            )
        },
    )
    exact_rhs = Matrix(
        [
            Rational(
                augmented_objective[equation].numerator,
                augmented_objective[equation].denominator,
            )
            for equation in selected_equations
        ]
    )
    solution_set = linsolve((exact_square, exact_rhs))
    solution_tuple = next(iter(solution_set))
    weights = [
        F(int(value.p), int(value.q)) for value in solution_tuple
    ]
    positive_inequality_weights = [
        (index, weight)
        for index, weight in enumerate(weights[:inequality_count])
        if weight > 0
    ]
    if positive_inequality_weights:
        largest_index, largest_weight = max(
            positive_inequality_weights, key=lambda item: item[1]
        )
        raise RuntimeError(
            "exact solution violates inequality dual signs: "
            f"{len(positive_inequality_weights)} positive, "
            f"largest row {largest_index} = {largest_weight} "
            f"(float {float(largest_weight)})"
        )

    exact_stationarity = [F(0)] * variable_count
    dual_objective = F(0)
    for weight, row, bound in zip(weights, rows, rhs, strict=True):
        for variable, coefficient in row.items():
            exact_stationarity[variable] += weight * coefficient
        dual_objective += weight * bound
    if exact_stationarity != objective:
        raise RuntimeError("full exact stationarity verification failed")
    if dual_objective != target_dual_objective:
        raise RuntimeError("full exact dual-objective verification failed")

    certificate = {
        "active_inequalities": [
            {
                "row": active["row"],
                "weight": str(weight),
            }
            for active, weight in zip(
                support["inequalities"],
                weights[:inequality_count],
                strict=True,
            )
            if weight
        ],
        "active_equalities": [
            {
                "row": active["row"],
                "weight": str(weight),
            }
            for active, weight in zip(
                support["equalities"],
                weights[inequality_count:],
                strict=True,
            )
            if weight
        ],
        "variable_objective": str(dual_objective),
    }
    payload = {
        "status": "perspective_global_anchor_basis_exact",
        "source": str(args.terminal_report),
        "anchor_search": str(args.anchor_search),
        "float_support": str(args.float_support),
        "diagnostics": diagnostics,
        "links": anchor["links"],
        "global_variable_objective": str(dual_objective),
        "exact_global_dual_certificate": certificate,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "global_variable_objective": str(dual_objective),
                "inequality_rows": len(certificate["active_inequalities"]),
                "equality_rows": len(certificate["active_equalities"]),
                "maximum_weight_denominator": max(
                    weight.denominator for weight in weights
                ),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
