#!/usr/bin/env python3
"""Expand exact games to full player orbits and solve all induced comparisons."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_intrinsic_permutation_union import permute_game
from n5_minkowski_family_product_search import component_margin, load_games


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source_games = load_games(args.input)
    permutations = list(itertools.permutations(range(5)))
    games = list(
        dict.fromkeys(
            permute_game(game, permutation)
            for game in source_games
            for permutation in permutations
        )
    )
    edges = set()
    add_induced_edges_fast(games, edges)
    best, orbit_games, quotient_edges, components = component_margin(games, edges)
    margin = None if best is None else best[0]
    found = margin is not None and margin < -1e-8
    output = {
        "status": (
            "equivariant_obstruction_found"
            if found
            else "no_equivariant_obstruction_found"
        ),
        "source": str(args.input),
        "margin": margin,
        "game_count": len(games),
        "edge_count": len(edges),
        "game_orbit_count": orbit_games,
        "quotient_edge_count": quotient_edges,
        "component_count": components,
        "family": {
            "games": [list(map(str, game)) for game in games],
            "edges": [list(edge) for edge in sorted(edges)],
        },
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: value for key, value in output.items() if key != "family"}, indent=2))
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
