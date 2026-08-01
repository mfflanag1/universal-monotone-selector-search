#!/usr/bin/env python3
"""Prune an equivariant family to games used by an optimal LP dual."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from n5_equivariant_quotient_lp import quotient_family, solve_quotient
from n5_minkowski_family_product_search import load_games


def components(game_count, edges):
    parent = list(range(game_count))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for lower, _lower_player, upper, _upper_player in edges:
        left = find(lower)
        right = find(upper)
        if left != right:
            parent[right] = left
    groups = {}
    for node in range(game_count):
        groups.setdefault(find(node), []).append(node)
    return sorted(groups.values(), key=len, reverse=True)


def restrict_family(games, edges, stabilizers, nodes):
    nodes = sorted(nodes)
    node_set = set(nodes)
    reindex = {node: index for index, node in enumerate(nodes)}
    return (
        [games[node] for node in nodes],
        {
            (reindex[lower], lower_player, reindex[upper], upper_player)
            for lower, lower_player, upper, upper_player in edges
            if lower in node_set and upper in node_set
        },
        [stabilizers[node] for node in nodes],
    )


def dual_games(result, game_count, edges, tolerance):
    core_row_count = 30 * game_count
    marginals = np.asarray(result.ineqlin.marginals)
    active = {
        row // 30
        for row in np.flatnonzero(np.abs(marginals[:core_row_count]) > tolerance)
    }
    ordered_edges = sorted(edges)
    for edge_index in np.flatnonzero(
        np.abs(marginals[core_row_count:]) > tolerance
    ):
        lower, _lower_player, upper, _upper_player = ordered_edges[int(edge_index)]
        active.update((lower, upper))
    return active


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    family = payload.get("family")
    if family is None:
        raise ValueError("input has no family")
    raw_games = load_games(args.input)
    raw_edges = {tuple(edge) for edge in family.get("edges", [])}
    games, edges, stabilizers = quotient_family(raw_games, raw_edges, 5)

    best = None
    for nodes in components(len(games), edges):
        local_games, local_edges, local_stabilizers = restrict_family(
            games, edges, stabilizers, nodes
        )
        margin = solve_quotient(local_games, local_edges, local_stabilizers, 5)
        if best is None or margin < best[0]:
            best = (margin, local_games, local_edges, local_stabilizers)
    if best is None:
        raise ValueError("input has no comparison component")

    margin, games, edges, stabilizers = best
    rounds = []
    while True:
        margin, result = solve_quotient(
            games, edges, stabilizers, 5, return_result=True
        )
        active = dual_games(
            result, len(games), edges, args.tolerance
        )
        row = {
            "games": len(games),
            "edges": len(edges),
            "dual_games": len(active),
            "margin": margin,
        }
        rounds.append(row)
        print(json.dumps(row), flush=True)
        if not active or len(active) == len(games):
            break
        candidate = restrict_family(games, edges, stabilizers, active)
        candidate_margin = solve_quotient(*candidate, 5)
        if abs(candidate_margin - margin) > 1e-8 * max(1.0, abs(margin)):
            row["rejected_candidate_margin"] = candidate_margin
            break
        games, edges, stabilizers = candidate

    output = {
        "status": "dual_support_pruned",
        "source": str(args.input),
        "margin": margin,
        "rounds": rounds,
        "family": {
            "games": [list(map(str, game)) for game in games],
            "quotient_edges": [list(edge) for edge in sorted(edges)],
        },
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
