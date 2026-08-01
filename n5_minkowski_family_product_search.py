#!/usr/bin/env python3
"""Search induced comparisons among sums or meets of two exact families."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_equivariant_quotient_lp import quotient_family, solve_quotient
from n5_intrinsic_permutation_union import permute_game


def load_games(path):
    payload = json.loads(path.read_text())
    family = payload.get("family")
    if family is None and isinstance(payload.get("best"), dict):
        family = payload["best"].get("family")
    if family is None:
        raise ValueError(f"no family in {path}")
    return [
        tuple(Fraction(value) for value in game)
        for game in family["games"]
    ]


def minimum_games(left, right):
    tables = [
        sorted({game[coalition] for game in left + right})
        for coalition in range(32)
    ]
    ranks = [
        {value: rank for rank, value in enumerate(table)}
        for table in tables
    ]
    left_codes = np.asarray(
        [
            [ranks[coalition][value] for coalition, value in enumerate(game)]
            for game in left
        ],
        dtype=np.uint16,
    )
    right_codes = np.asarray(
        [
            [ranks[coalition][value] for coalition, value in enumerate(game)]
            for game in right
        ],
        dtype=np.uint16,
    )
    combined = np.minimum(
        left_codes[:, np.newaxis, :], right_codes[np.newaxis, :, :]
    ).reshape(-1, 32)
    unique = np.unique(combined, axis=0)
    return [
        tuple(tables[coalition][int(code)] for coalition, code in enumerate(row))
        for row in unique
    ]


def component_margin(games, edges):
    canonical_games, quotient_edges, stabilizers = quotient_family(
        games, edges, 5
    )
    parent = list(range(len(canonical_games)))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for lower, _lower_player, upper, _upper_player in quotient_edges:
        left = find(lower)
        right = find(upper)
        if left != right:
            parent[right] = left
    components = {}
    for node in range(len(canonical_games)):
        components.setdefault(find(node), []).append(node)
    best = None
    for component_index, nodes in enumerate(
        sorted(components.values(), key=len, reverse=True), start=1
    ):
        node_set = set(nodes)
        reindex = {node: index for index, node in enumerate(nodes)}
        local_edges = {
            (reindex[lower], lower_player, reindex[upper], upper_player)
            for lower, lower_player, upper, upper_player in quotient_edges
            if lower in node_set and upper in node_set
        }
        local_games = [canonical_games[node] for node in nodes]
        local_stabilizers = [stabilizers[node] for node in nodes]
        margin = solve_quotient(
            local_games, local_edges, local_stabilizers, 5
        )
        if best is None or margin < best[0]:
            best = (margin, local_games, local_edges)
            print(
                json.dumps(
                    {
                        "component": component_index,
                        "components": len(components),
                        "games": len(local_games),
                        "edges": len(local_edges),
                        "margin": margin,
                    }
                ),
                flush=True,
            )
        if margin < -1e-8:
            break
    return best, len(canonical_games), len(quotient_edges), len(components)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", type=Path, required=True)
    parser.add_argument("--right", type=Path)
    parser.add_argument("--permutation", default="0,1,2,3,4")
    parser.add_argument("--operation", choices=("sum", "minimum"), default="sum")
    parser.add_argument("--save-generated", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    permutation = tuple(int(value) for value in args.permutation.split(","))
    if sorted(permutation) != list(range(5)):
        parser.error("--permutation must list 0,1,2,3,4 exactly once")
    left = load_games(args.left)
    right = load_games(args.right or args.left)
    right = [permute_game(game, permutation) for game in right]
    left_grands = {game[-1] for game in left}
    right_grands = {game[-1] for game in right}
    if len(left_grands) != 1 or len(right_grands) != 1:
        raise ValueError("each input family must have a common grand worth")
    if args.operation == "minimum":
        left_grand = next(iter(left_grands))
        right_grand = next(iter(right_grands))
        if not right_grand:
            raise ValueError("the right grand worth must be nonzero")
        scale = left_grand / right_grand
        right = [
            tuple(scale * value for value in game) for game in right
        ]
    if args.operation == "minimum":
        games = minimum_games(left, right)
    else:
        games = list(
            dict.fromkeys(
                tuple(
                    left_value + right_value
                    for left_value, right_value in zip(a, b)
                )
                for a in left
                for b in right
            )
        )
    edges = set()
    add_induced_edges_fast(games, edges)
    print(
        json.dumps(
            {
                "phase": "induced",
                "games": len(games),
                "edges": len(edges),
            }
        ),
        flush=True,
    )
    best, orbit_games, quotient_edges, component_count = component_margin(
        games, edges
    )
    margin = None if best is None else best[0]
    found = margin is not None and margin < -1e-8
    payload = {
        "status": (
            "equivariant_obstruction_found"
            if found
            else "no_equivariant_obstruction_found"
        ),
        "left": str(args.left),
        "right": str(args.right or args.left),
        "permutation": list(permutation),
        "operation": args.operation,
        "sum_game_count": len(games),
        "induced_edge_count": len(edges),
        "game_orbit_count": orbit_games,
        "quotient_edge_count": quotient_edges,
        "component_count": component_count,
        "margin": margin,
        "quotient_component_games": (
            [list(map(str, game)) for game in best[1]] if found else None
        ),
        "quotient_component_edges": (
            [list(edge) for edge in sorted(best[2])] if found else None
        ),
    }
    if args.save_generated:
        payload["family"] = {
            "games": [list(map(str, game)) for game in games],
            "edges": [list(edge) for edge in sorted(edges)],
        }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if not key.startswith("quotient_component_") and key != "family"}, indent=2))
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
