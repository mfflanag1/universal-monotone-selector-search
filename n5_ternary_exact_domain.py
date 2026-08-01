#!/usr/bin/env python3
"""Enumerate the complete {0,1,2}-valued monotone exact n=5 domain."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_boolean_exact_domain import GRAND, N, antichains, game_code
from n5_facet_search import load_facets
from n5_point_deletion_family_search import connected_components
from n5_sparse_family_lp import sparse_family_slack_float


def monotone_boolean_codes() -> list[int]:
    return [
        game_code(antichain)
        for antichain in antichains()
        if antichain
    ]


def vector(code: int) -> np.ndarray:
    return np.asarray([(code >> coalition) & 1 for coalition in range(GRAND + 1)], dtype=np.int8)


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facets, _types = load_facets()
    facet_matrix = np.asarray(facets, dtype=np.int16)[:, 1:]
    codes = monotone_boolean_codes()
    boolean_games = np.stack([vector(code) for code in codes])
    scores = boolean_games[:, 1:].astype(np.int16) @ facet_matrix.T
    exact_games: set[bytes] = {bytes(GRAND + 1)}
    exact_boolean = np.all(scores >= 0, axis=1)
    for index in np.flatnonzero(exact_boolean):
        exact_games.add(bytes(boolean_games[index]))

    for upper_index, upper_code in enumerate(codes):
        nested = np.fromiter(
            (index for index, lower_code in enumerate(codes) if lower_code & ~upper_code == 0),
            dtype=np.int32,
        )
        feasible = nested[
            np.all(scores[nested] + scores[upper_index] >= 0, axis=1)
        ]
        upper = boolean_games[upper_index]
        for lower_index in feasible:
            exact_games.add(bytes(upper + boolean_games[lower_index]))
        if upper_index % 500 == 0:
            print(json.dumps({"phase": "games", "upper": upper_index, "exact": len(exact_games)}), flush=True)

    games = sorted(exact_games)
    index = {game: position for position, game in enumerate(games)}
    edges = set()
    for lower, game in enumerate(games):
        for coalition in range(1, GRAND + 1):
            value = game[coalition]
            for raised in range(value + 1, 3):
                successor = bytearray(game)
                successor[coalition] = raised
                upper = index.get(bytes(successor))
                if upper is not None:
                    edges.add((lower, upper, coalition))
    components = connected_components(len(games), edges)
    rows = []
    best = None
    for component, (nodes, local_edges) in enumerate(components):
        local_games = [tuple(games[node]) for node in nodes]
        margin = sparse_family_slack_float(local_games, local_edges, N, "highs-ipm")
        row = {
            "component": component,
            "nodes": len(nodes),
            "edges": len(local_edges),
            "cycle_rank": len(local_edges) - len(nodes) + 1,
            "margin": margin,
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
        if margin is not None and (best is None or margin < best[0]):
            best = (margin, local_games, local_edges)
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "scope": "complete monotone {0,1,2}-valued five-player exact-game domain",
        "monotone_boolean_threshold_count": len(codes),
        "exact_game_count": len(games),
        "legal_edge_count": len(edges),
        "component_count": len(components),
        "components": rows,
        "best": None
        if best is None
        else {
            "margin": best[0],
            "games": [list(game) for game in best[1]] if best[0] < -1e-8 else None,
            "edges": [list(edge) for edge in best[2]] if best[0] < -1e-8 else None,
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if best is not None and best[0] < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
