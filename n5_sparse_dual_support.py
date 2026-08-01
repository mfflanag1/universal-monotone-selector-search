#!/usr/bin/env python3
"""Extract an induced family from a large sparse LP's active dual rows."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from family_search import (
    Family,
    common_core_gap_float,
    members,
    serial_family,
)
from n5_m5_bridge_closure import induced_edges
from n5_sparse_family_lp import (
    sparse_family_result,
    sparse_family_slack_float,
)


F = Fraction


def dual_envelope_support(
    all_games,
    all_edges,
    n: int,
    result,
    source: str,
    threshold: float,
    method: str,
):
    grand = (1 << n) - 1
    active_nodes = set()
    active_core_rows = 0
    active_monotonicity_rows = 0
    row = 0
    marginals = result.ineqlin.marginals
    for game_index in range(len(all_games)):
        for _coalition in range(1, grand):
            if abs(float(marginals[row])) > threshold:
                active_nodes.add(game_index)
                active_core_rows += 1
            row += 1
    for lower, upper, coalition in all_edges:
        for _player in members(coalition, n):
            if abs(float(marginals[row])) > threshold:
                active_nodes.update((lower, upper))
                active_monotonicity_rows += 1
            row += 1
    if row != len(marginals):
        raise RuntimeError("dual row count mismatch")

    envelope_nodes = set()
    for coalition in range(1, grand):
        maximum = max(game[coalition] for game in all_games)
        for node, game in enumerate(all_games):
            if game[coalition] == maximum:
                envelope_nodes.add(node)
                break
    source_nodes = sorted(active_nodes | envelope_nodes)
    games = [all_games[node] for node in source_nodes]
    edges = induced_edges(games, grand)
    if any(game[grand] != games[0][grand] for game in games):
        grand_groups = defaultdict(list)
        for node, game in enumerate(games):
            grand_groups[game[:grand]].append(node)
        aggregate_edges = set()
        for nodes in grand_groups.values():
            ordered = sorted(
                nodes, key=lambda node: games[node][grand]
            )
            for position, lower in enumerate(ordered):
                for upper in ordered[position + 1 :]:
                    if games[lower][grand] < games[upper][grand]:
                        aggregate_edges.add((lower, upper, grand))
        edges = sorted(set(edges) | aggregate_edges)
    family = Family(
        games,
        edges,
        [0] * len(games),
        "n5_sparse_dual_envelope_support",
    )
    margin = sparse_family_slack_float(
        games, edges, n, method=method
    )
    common_gap = (
        common_core_gap_float(games, n)
        if all(
            game[grand] == games[0][grand] for game in games
        )
        else None
    )
    return {
        "status": "sparse_dual_envelope_support_complete",
        "n": n,
        "source": source,
        "source_margin_float": float(result.x[len(all_games) * n]),
        "active_core_rows": active_core_rows,
        "active_monotonicity_rows": active_monotonicity_rows,
        "active_node_count": len(active_nodes),
        "envelope_node_count": len(envelope_nodes),
        "source_nodes": source_nodes,
        "node_count": len(games),
        "edge_count": len(edges),
        "common_core_budget_gap_float": common_gap,
        "max_min_monotonicity_margin_float": margin,
        "family": serial_family(family, margin, common_gap),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=1e-9)
    parser.add_argument(
        "--method",
        choices=("highs", "highs-ds", "highs-ipm"),
        default="highs-ipm",
    )
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    n = int(payload["n"])
    raw = payload["family"]
    all_games = [
        tuple(F(value) for value in game)
        for game in raw["games"]
    ]
    all_edges = [
        (
            int(edge["lower"]),
            int(edge["upper"]),
            int(edge["coalition"]),
        )
        for edge in raw["edges"]
    ]
    result = sparse_family_result(
        all_games, all_edges, n, args.method
    )
    if not result.success:
        raise RuntimeError(result.message)

    result_payload = dual_envelope_support(
        all_games,
        all_edges,
        n,
        result,
        str(args.input),
        args.threshold,
        args.method,
    )
    args.output.write_text(
        json.dumps(result_payload, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                key: result_payload[key]
                for key in (
                    "status",
                    "source_margin_float",
                    "active_core_rows",
                    "active_monotonicity_rows",
                    "active_node_count",
                    "envelope_node_count",
                    "node_count",
                    "edge_count",
                    "common_core_budget_gap_float",
                    "max_min_monotonicity_margin_float",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
