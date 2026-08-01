#!/usr/bin/env python3
"""Close an exact family under meets with a small set of generators."""

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
    parser.add_argument("--generators", type=int, default=5)
    parser.add_argument("--indices")
    parser.add_argument("--solve-every", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--save-family", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    permutation = tuple(int(value) for value in args.permutation.split(","))
    if sorted(permutation) != list(range(5)):
        parser.error("--permutation must list 0,1,2,3,4 exactly once")
    games = load_games(args.left)
    right = [
        permute_game(game, permutation) for game in load_games(args.right)
    ]
    left_grands = {game[-1] for game in games}
    right_grands = {game[-1] for game in right}
    if len(left_grands) != 1 or len(right_grands) != 1:
        raise ValueError("each input family must have one grand worth")
    scale = next(iter(left_grands)) / next(iter(right_grands))
    right = [tuple(scale * value for value in game) for game in right]
    if args.indices:
        indices = [int(value) for value in args.indices.split(",")]
        right = [right[index] for index in indices]
    else:
        random.Random(args.seed).shuffle(right)

    game_set = set(games)
    rows = []
    best = None
    edges = set()
    for generator_index, generator in enumerate(
        right[: args.generators], start=1
    ):
        prior_games = list(games)
        for game in prior_games:
            meet = tuple(
                min(left_value, right_value)
                for left_value, right_value in zip(game, generator)
            )
            if meet not in game_set:
                game_set.add(meet)
                games.append(meet)
        growth = len(games) - len(prior_games)
        print(
            json.dumps(
                {
                    "phase": "generator",
                    "generator": generator_index,
                    "games": len(games),
                    "growth": growth,
                }
            ),
            flush=True,
        )
        if generator_index % args.solve_every and generator_index != args.generators:
            continue
        edges.clear()
        add_induced_edges_fast(games, edges)
        component, orbit_games, quotient_edges, components = component_margin(
            games, edges
        )
        margin = None if component is None else component[0]
        row = {
            "generator": generator_index,
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
    output = {
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
    }
    if args.save_family or found:
        output["family"] = {
            "games": [list(map(str, game)) for game in games],
            "edges": [list(edge) for edge in sorted(edges)],
        }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
