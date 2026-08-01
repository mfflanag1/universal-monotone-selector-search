#!/usr/bin/env python3
"""Search the symmetry quotient of a growing random point-generated exact domain."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from pathlib import Path

from mixed_point_family_fast import random_points
from n5_bppv_pool_union_search import games_from_sample
from n5_equivariant_quotient_lp import equivariant_quotient_margin
from n5_point_deletion_walk_search import add_induced_edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batches", type=int, default=100)
    parser.add_argument("--pool-count", type=int, default=24)
    parser.add_argument("--sample-count", type=int, default=9)
    parser.add_argument("--total", type=int, default=9)
    parser.add_argument("--bound", type=int, default=2)
    parser.add_argument(
        "--pool-kind",
        choices=("random", "bppv", "housman-clark", "mixed"),
        default="random",
    )
    parser.add_argument("--anchor", action="store_true")
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--solve-every", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    n = 5
    grand = (1 << n) - 1
    bppv_seed = ((3, 1, 1, 2, 2), (2, 3, 2, 1, 1))
    housman_clark_seed = (
        (0, 0, 4, 4, 0),
        (0, 1, 4, 3, 0),
        (1, 0, 3, 4, 0),
        (3, 4, 1, 0, 0),
        (4, 3, 0, 1, 0),
        (4, 4, 0, 0, 0),
    )
    anchors = []
    if args.pool_kind == "random":
        points = random_points(rng, n, args.pool_count, args.total, args.bound)
    else:
        if args.pool_kind == "bppv":
            seed = bppv_seed
        elif args.pool_kind == "housman-clark":
            seed = housman_clark_seed
        else:
            seed = bppv_seed + tuple(
                point[:-1] + (point[-1] + 1,)
                for point in housman_clark_seed
            )
        points = sorted(
            {
                tuple(point[index] for index in permutation)
                for point in seed
                for permutation in itertools.permutations(range(n))
            }
        )
        if args.anchor:
            anchors = list(seed)
    if len(anchors) > args.sample_count:
        parser.error("the structured anchors exceed --sample-count")
    point_values = [
        tuple(
            sum(
                point[player]
                for player in range(n)
                if coalition >> player & 1
            )
            for coalition in range(grand + 1)
        )
        for point in points
    ]
    games = []
    game_set = set()
    edges = set()
    margin = None
    orbit_games = 0
    quotient_edges = 0
    for batch in range(1, args.batches + 1):
        anchor_indices = [points.index(point) for point in anchors]
        choices = [index for index in range(len(points)) if index not in anchor_indices]
        sample = anchor_indices + rng.sample(
            choices, args.sample_count - len(anchor_indices)
        )
        for game in games_from_sample(point_values, sample):
            if game not in game_set:
                game_set.add(game)
                games.append(game)
        if batch % args.solve_every and batch != args.batches:
            continue
        edges.clear()
        add_induced_edges(games, edges)
        margin, orbit_games, quotient_edges = equivariant_quotient_margin(
            games, edges, n
        )
        print(
            json.dumps(
                {
                    "batch": batch,
                    "games": len(games),
                    "edges": len(edges),
                    "game_orbits": orbit_games,
                    "quotient_edges": quotient_edges,
                    "margin": margin,
                }
            ),
            flush=True,
        )
        if margin < -1e-8:
            break
    found = margin is not None and margin < -1e-8
    payload = {
        "status": (
            "equivariant_obstruction_found"
            if found
            else "no_equivariant_obstruction_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "points": [list(point) for point in points],
        "margin": margin,
        "game_orbits": orbit_games,
        "quotient_edges": quotient_edges,
        "games": [list(game) for game in games] if found else None,
        "edges": [list(edge) for edge in sorted(edges)] if found else None,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"status": payload["status"], "margin": margin}))
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
