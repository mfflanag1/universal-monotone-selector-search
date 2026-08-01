#!/usr/bin/env python3
"""Mix an exact family with a shortest-path slice of another exact family."""

from __future__ import annotations

import argparse
import collections
import json
from fractions import Fraction
from pathlib import Path

from family_search import common_core_gap_float
from n5_facet_search import margin_dual, permute_mask


F = Fraction
N = 5
GRAND = 31


def load_archive(path: Path) -> tuple[list[list[F]], list[tuple[int, int, int, F]]]:
    family = json.loads(path.read_text())["family"]
    games = [[F(value) for value in game] for game in family["games"]]
    edges = [
        (
            int(edge["lower"]),
            int(edge["upper"]),
            int(edge["coalition"]),
            F(edge["delta"]),
        )
        for edge in family["edges"]
    ]
    return games, edges


def shortest_path(
    node_count: int,
    edges: list[tuple[int, int, int, F]],
    sources: list[int],
    targets: set[int],
) -> list[int]:
    adjacency = [[] for _ in range(node_count)]
    for lower, upper, _, _ in edges:
        adjacency[lower].append(upper)
        adjacency[upper].append(lower)
    queue = collections.deque(sources)
    predecessor = {source: None for source in sources}
    found = None
    while queue:
        node = queue.popleft()
        if node in targets:
            found = node
            break
        for neighbor in adjacency[node]:
            if neighbor in predecessor:
                continue
            predecessor[neighbor] = node
            queue.append(neighbor)
    if found is None:
        raise ValueError("source and target sets are disconnected")
    path = []
    current = found
    while current is not None:
        path.append(current)
        current = predecessor[current]
    path.reverse()
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--mixer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mixer-sources", required=True)
    parser.add_argument("--mixer-targets", required=True)
    parser.add_argument("--permutation", required=True)
    parser.add_argument("--scale", type=F, required=True)
    args = parser.parse_args()

    base_games, base_edges = load_archive(args.base)
    mixer_games, mixer_edges = load_archive(args.mixer)
    sources = [int(value) for value in args.mixer_sources.split(",")]
    targets = {
        int(value) for value in args.mixer_targets.split(",")
    }
    permutation = tuple(
        int(value) for value in args.permutation.split(",")
    )
    if sorted(permutation) != list(range(N)):
        raise ValueError("permutation must contain 0,1,2,3,4")
    path = shortest_path(
        len(mixer_games), mixer_edges, sources, targets
    )
    path_set = set(path)
    mixer_node_map = {
        old_node: new_node for new_node, old_node in enumerate(path)
    }
    sliced_edges = [
        (
            mixer_node_map[lower],
            mixer_node_map[upper],
            permute_mask(coalition, permutation),
            delta * args.scale,
        )
        for lower, upper, coalition, delta in mixer_edges
        if lower in path_set and upper in path_set
    ]
    sliced_games = []
    for old_node in path:
        game = [F(0)] * (GRAND + 1)
        for coalition in range(GRAND + 1):
            game[permute_mask(coalition, permutation)] = (
                mixer_games[old_node][coalition] * args.scale
            )
        sliced_games.append(game)

    product_games = []
    for base_game in base_games:
        for mixer_game in sliced_games:
            product_games.append(
                [
                    base_game[coalition] + mixer_game[coalition]
                    for coalition in range(GRAND + 1)
                ]
            )
    layer_size = len(sliced_games)

    product_edges = []
    product_bumps = []
    for mixer_node in range(layer_size):
        for lower, upper, coalition, delta in base_edges:
            product_edges.append(
                (
                    lower * layer_size + mixer_node,
                    upper * layer_size + mixer_node,
                    coalition,
                )
            )
            product_bumps.append(float(delta))
    for base_node in range(len(base_games)):
        for lower, upper, coalition, delta in sliced_edges:
            product_edges.append(
                (
                    base_node * layer_size + lower,
                    base_node * layer_size + upper,
                    coalition,
                )
            )
            product_bumps.append(float(delta))

    games_float = [
        [float(value) for value in game] for game in product_games
    ]
    margin, _, _, _ = margin_dual(games_float, product_edges)
    gap = common_core_gap_float(games_float, N)
    result = {
        "status": "mixed_product_path_float",
        "base": str(args.base),
        "mixer": str(args.mixer),
        "mixer_path": path,
        "mixer_induced_edge_count": len(sliced_edges),
        "permutation": list(permutation),
        "scale": str(args.scale),
        "common_core_gap_float": gap,
        "record": {
            "n": N,
            "edges": [list(edge) for edge in product_edges],
            "bumps": product_bumps,
            "best_margin_float": margin,
            "best_games_float": games_float,
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "mixer_path": path,
                "mixer_path_nodes": len(path),
                "mixer_induced_edges": len(sliced_edges),
                "games": len(product_games),
                "edges": len(product_edges),
                "margin": margin,
                "common_core_gap": gap,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
