#!/usr/bin/env python3
"""Extract the protected-flow dual from a BPPV pool-union family."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from n5_mixed_terminal_search import lower_expectation_dual
from n5_point_deletion_family_search import connected_components
from n5_sparse_family_lp import sparse_family_result


F = Fraction
N = 5
GRAND = (1 << N) - 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--component-limit", type=int, default=20)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    games = [tuple(F(value) for value in game) for game in payload["games"]]
    edges = {tuple(int(value) for value in edge) for edge in payload["edges"]}
    components = connected_components(len(games), edges)
    best = None
    rows = []
    for component_index, (nodes, local_edges) in enumerate(
        components[: args.component_limit]
    ):
        local_games = [games[node] for node in nodes]
        result = sparse_family_result(
            local_games, local_edges, N, method="highs-ipm"
        )
        if not result.success:
            continue
        margin = float(result.x[len(local_games) * N])
        rows.append(
            {
                "component": component_index,
                "nodes": len(nodes),
                "edges": len(local_edges),
                "margin": margin,
            }
        )
        if best is None or margin < best[0]:
            best = (margin, component_index, nodes, local_edges, result)
    if best is None:
        raise RuntimeError("no component LP solved")
    margin, component_index, nodes, local_edges, result = best
    local_games = [games[node] for node in nodes]
    core_row_count = len(local_games) * (GRAND - 1)
    marginals = result.ineqlin.marginals[core_row_count:]
    divergence = [[0.0] * N for _ in local_games]
    active_edges = []
    row = 0
    for edge_index, (lower, upper, coalition) in enumerate(local_edges):
        for player in range(N):
            if not coalition >> player & 1:
                continue
            weight = -float(marginals[row])
            row += 1
            if weight <= 1e-9:
                continue
            divergence[lower][player] += weight
            divergence[upper][player] -= weight
            active_edges.append(
                {
                    "edge": edge_index,
                    "lower": lower,
                    "upper": upper,
                    "coalition": coalition,
                    "player": player,
                    "weight": weight,
                }
            )
    active_nodes = [index for index, vector in enumerate(divergence) if any(abs(value) > 1e-9 for value in vector)]
    lower_rows = []
    for node in active_nodes:
        dual = lower_expectation_dual(
            local_games[node], tuple(divergence[node]), N
        )
        if dual is None:
            inequalities = []
            rhs = []
            for coalition in range(1, GRAND):
                inequalities.append(
                    [
                        -1.0 if coalition >> player & 1 else 0.0
                        for player in range(N)
                    ]
                )
                rhs.append(-float(local_games[node][coalition]))
            primal = linprog(
                divergence[node],
                A_ub=np.asarray(inequalities),
                b_ub=np.asarray(rhs),
                A_eq=np.ones((1, N)),
                b_eq=np.asarray([float(local_games[node][GRAND])]),
                bounds=[(None, None)] * N,
                method="highs",
            )
            if not primal.success:
                raise RuntimeError("lower-expectation primal and dual failed")
            value = float(primal.fun)
            coalition_weights = []
            beta = None
        else:
            value, coalition_weights, beta = dual
        lower_rows.append(
            {
                "node": node,
                "original_node": nodes[node],
                "divergence": divergence[node],
                "lower_expectation": value,
                "beta": beta,
                "coalition_weights": [
                    {"coalition": coalition, "weight": float(weight)}
                    for coalition, weight in enumerate(coalition_weights, start=1)
                    if weight > 1e-9
                ],
                "game": [str(value) for value in local_games[node]],
            }
        )
    report = {
        "status": "bppv_pool_active_dual_report",
        "source": str(args.input),
        "component_rows": rows,
        "selected_component": component_index,
        "node_count": len(nodes),
        "edge_count": len(local_edges),
        "margin": margin,
        "active_monotonicity_rows": len(active_edges),
        "active_nodes": len(active_nodes),
        "monotonicity_weight_sum": sum(row["weight"] for row in active_edges),
        "lower_expectation_sum": sum(row["lower_expectation"] for row in lower_rows),
        "active_edges": active_edges,
        "lower_expectations": lower_rows,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {key: report[key] for key in report if key not in ("active_edges", "lower_expectations", "component_rows")},
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
