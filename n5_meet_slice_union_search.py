#!/usr/bin/env python3
"""Grow an exact meet family one complete right-factor slice at a time."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_intrinsic_permutation_union import permute_game
from n5_minkowski_family_product_search import component_margin, load_games


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", type=Path, required=True)
    parser.add_argument("--right", type=Path, required=True)
    parser.add_argument("--permutation", default="0,1,2,3,4")
    parser.add_argument("--batches", type=int, default=30)
    parser.add_argument("--solve-every", type=int, default=1)
    parser.add_argument("--include-left", action="store_true")
    parser.add_argument("--save-family", action="store_true")
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    permutation = tuple(int(value) for value in args.permutation.split(","))
    left = load_games(args.left)
    right = [
        permute_game(game, permutation) for game in load_games(args.right)
    ]
    left_grand = {game[-1] for game in left}
    right_grand = {game[-1] for game in right}
    if len(left_grand) != 1 or len(right_grand) != 1:
        raise ValueError("each input must have one grand worth")
    scale = next(iter(left_grand)) / next(iter(right_grand))
    right = [tuple(scale * value for value in game) for game in right]
    rng = random.Random(args.seed)
    rng.shuffle(right)
    games = list(left) if args.include_left else []
    game_set = set(games)
    rows = []
    best = None
    edges = set()
    for batch, right_game in enumerate(right[: args.batches], start=1):
        for left_game in left:
            game = tuple(
                min(left_value, right_value)
                for left_value, right_value in zip(left_game, right_game)
            )
            if game not in game_set:
                game_set.add(game)
                games.append(game)
        if batch % args.solve_every and batch != args.batches:
            continue
        edges.clear()
        add_induced_edges_fast(games, edges)
        component, orbit_games, quotient_edges, components = component_margin(
            games, edges
        )
        margin = None if component is None else component[0]
        row = {
            "batch": batch,
            "games": len(games),
            "edges": len(edges),
            "game_orbits": orbit_games,
            "quotient_edges": quotient_edges,
            "components": components,
            "margin": margin,
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
        if best is None or (margin is not None and margin < best[0]):
            best = (margin, component)
        if margin is not None and margin < -1e-8:
            break
    found = best is not None and best[0] is not None and best[0] < -1e-8
    payload = {
        "status": (
            "equivariant_obstruction_found"
            if found
            else "no_equivariant_obstruction_found"
        ),
        "configuration": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
        },
        "rows": rows,
        "margin": None if best is None else best[0],
        "games": [list(map(str, game)) for game in games] if found else None,
        "edges": [list(edge) for edge in sorted(edges)] if found else None,
    }
    if args.save_family:
        payload["family"] = {
            "games": [list(map(str, game)) for game in games],
            "edges": [list(edge) for edge in sorted(edges)],
        }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
