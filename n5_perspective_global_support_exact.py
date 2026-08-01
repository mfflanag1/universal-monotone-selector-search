#!/usr/bin/env python3
"""Exactify the global block of a perspective dual support."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from n5_facet_search import load_facets
from n5_terminal_topology_relaxation import reconstruct_exact_dual


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("float_support", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    support = json.loads(args.float_support.read_text())
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    node_index = {node: index for index, node in enumerate(node_ids)}
    base_variable_count = len(node_ids) * coalition_count
    mixed_nodes = {
        int(terminal["node"])
        for terminal in terminals
        if terminal.get("coefficient") is None
        and len(set(terminal["scaled_divergence"])) > 2
    }

    def column(node: int, coalition: int) -> int:
        return node_index[node] * coalition_count + coalition

    protected_players: dict[tuple[int, int], set[int]] = defaultdict(set)
    for path in report["terminal_path_decomposition"]:
        protected_players[
            int(path["source"]), int(path["sink"])
        ].add(int(path["player"]))
    facets, _ = load_facets()

    def row_from_metadata(metadata: list[Any]) -> tuple[dict[int, F], F]:
        row_type = metadata[0]
        if row_type == "exact_cone":
            _, node, facet_index = metadata
            facet = facets[int(facet_index)]
            return (
                {
                    column(int(node), coalition): -F(coefficient)
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                },
                F(0),
            )
        if row_type == "directed_monotonicity":
            _, source, sink, coalition = metadata
            return (
                {
                    column(int(source), int(coalition)): F(1),
                    column(int(sink), int(coalition)): F(-1),
                },
                F(0),
            )
        if row_type == "empty":
            _, node = metadata
            return ({column(int(node), 0): F(1)}, F(0))
        if row_type == "grand":
            _, node = metadata
            return ({column(int(node), grand): F(1)}, F(1))
        if row_type == "protected_invariance":
            _, source, sink, coalition = metadata
            return (
                {
                    column(int(source), int(coalition)): F(1),
                    column(int(sink), int(coalition)): F(-1),
                },
                F(0),
            )
        raise RuntimeError(f"unsupported global support row {row_type}")

    inequality_entries = [
        active
        for active in support["inequalities"]
        if not active["row"][0].startswith("perspective_")
    ]
    equality_entries = [
        active
        for active in support["equalities"]
        if not active["row"][0].startswith("perspective_")
    ]
    inequality_rows_full = [
        row_from_metadata(active["row"])[0]
        for active in inequality_entries
    ]
    equality_rows_full = [
        row_from_metadata(active["row"])[0]
        for active in equality_entries
    ]
    inequality_rhs = [
        row_from_metadata(active["row"])[1]
        for active in inequality_entries
    ]
    equality_rhs = [
        row_from_metadata(active["row"])[1]
        for active in equality_entries
    ]

    kept_full_columns = [
        column(node, coalition)
        for node in node_ids
        if node not in mixed_nodes
        for coalition in range(coalition_count)
    ]
    compressed = {
        full_column: index
        for index, full_column in enumerate(kept_full_columns)
    }
    inequality_rows = [
        {
            compressed[full_column]: float(value)
            for full_column, value in row.items()
            if full_column in compressed
        }
        for row in inequality_rows_full
    ]
    equality_rows = [
        {
            compressed[full_column]: float(value)
            for full_column, value in row.items()
            if full_column in compressed
        }
        for row in equality_rows_full
    ]

    objective_full = [F(0)] * base_variable_count
    for terminal in terminals:
        node = int(terminal["node"])
        if node in mixed_nodes:
            continue
        raw_coefficient = terminal.get("coefficient")
        if raw_coefficient is not None:
            coefficient = F(raw_coefficient)
            coalition = int(terminal["coalition"])
            if coefficient > 0:
                objective_full[column(node, coalition)] -= coefficient
            elif coefficient < 0:
                objective_full[
                    column(node, grand ^ coalition)
                ] -= -coefficient
        else:
            divergence = [F(value) for value in terminal["scaled_divergence"]]
            low, high = sorted(set(divergence))
            coalition = sum(
                1 << player
                for player, value in enumerate(divergence)
                if value == high
            )
            objective_full[column(node, coalition)] -= high - low
            objective_full[column(node, grand)] -= low
    objective = [objective_full[full] for full in kept_full_columns]

    support_rows = inequality_rows + equality_rows
    row_indices = []
    columns = []
    values = []
    for support_column, row in enumerate(support_rows):
        for variable, coefficient in row.items():
            row_indices.append(variable)
            columns.append(support_column)
            values.append(coefficient)
    matrix = coo_matrix(
        (values, (row_indices, columns)),
        shape=(len(objective), len(support_rows)),
    ).tocsr()
    bounds = (
        [(None, 0.0)] * len(inequality_rows)
        + [(None, None)] * len(equality_rows)
    )
    solve = linprog(
        np.zeros(len(support_rows)),
        A_eq=matrix,
        b_eq=np.asarray([float(value) for value in objective]),
        bounds=bounds,
        method="highs-ds",
    )
    if not solve.success:
        raise RuntimeError(solve.message)
    pseudo_result = SimpleNamespace(
        ineqlin=SimpleNamespace(
            marginals=np.asarray(solve.x[: len(inequality_rows)])
        ),
        eqlin=SimpleNamespace(
            marginals=np.asarray(solve.x[len(inequality_rows) :])
        ),
    )
    exact_dual = reconstruct_exact_dual(
        pseudo_result,
        inequality_rows,
        [float(value) for value in inequality_rhs],
        [tuple(active["row"]) for active in inequality_entries],
        equality_rows,
        [float(value) for value in equality_rhs],
        [tuple(active["row"]) for active in equality_entries],
        objective,
        None,
        active_tolerance=1e-12,
    )
    if exact_dual is None:
        raise RuntimeError("failed to exactify the projected global support")

    weights = {
        tuple(active["row"]): F(active["weight"])
        for active in exact_dual["active_inequalities"]
        + exact_dual["active_equalities"]
    }
    stationarity = [F(0)] * base_variable_count
    global_objective = F(0)
    for active, row, rhs in zip(
        inequality_entries,
        inequality_rows_full,
        inequality_rhs,
        strict=True,
    ):
        weight = weights.get(tuple(active["row"]), F(0))
        for variable, coefficient in row.items():
            stationarity[variable] += weight * coefficient
        global_objective += weight * rhs
    for active, row, rhs in zip(
        equality_entries,
        equality_rows_full,
        equality_rhs,
        strict=True,
    ):
        weight = weights.get(tuple(active["row"]), F(0))
        for variable, coefficient in row.items():
            stationarity[variable] += weight * coefficient
        global_objective += weight * rhs
    links = {
        node: [
            -stationarity[column(node, coalition)]
            for coalition in range(coalition_count)
        ]
        for node in mixed_nodes
    }
    if any(
        stationarity[full_column] != objective_full[full_column]
        for full_column in kept_full_columns
    ):
        raise RuntimeError("projected global stationarity failed")

    payload = {
        "status": "perspective_global_support_exact",
        "source": str(args.terminal_report),
        "float_support": str(args.float_support),
        "global_variable_objective": str(global_objective),
        "links": {
            str(node): [str(value) for value in values]
            for node, values in links.items()
        },
        "exact_projected_dual_certificate": exact_dual,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "global_variable_objective": payload[
                    "global_variable_objective"
                ],
                "link_denominator_max": max(
                    value.denominator
                    for values in links.values()
                    for value in values
                ),
                "dual_inequality_rows": len(
                    exact_dual["active_inequalities"]
                ),
                "dual_equality_rows": len(
                    exact_dual["active_equalities"]
                ),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
