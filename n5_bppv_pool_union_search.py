#!/usr/bin/env python3
"""Accumulate legal exact games from subsets of the full BPPV point orbit."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from pathlib import Path

from n5_point_deletion_family_search import connected_components
from n5_point_deletion_walk_search import (
    GRAND,
    N,
    add_induced_edges,
    orbit_points,
)
from n5_intrinsic_permutation_union import permute_game
from n5_sparse_family_lp import sparse_family_slack_float


def games_from_sample(point_values, sample):
    games = set()
    for state in range(1, 1 << len(sample)):
        active = [
            sample[index]
            for index in range(len(sample))
            if state >> index & 1
        ]
        games.add(
            tuple(
                min(point_values[index][coalition] for index in active)
                for coalition in range(GRAND + 1)
            )
        )
    return games


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--point-count", type=int, default=10)
    parser.add_argument("--permutation-close", action="store_true")
    parser.add_argument(
        "--method", choices=("highs", "highs-ds", "highs-ipm"),
        default="highs-ipm"
    )
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    points = orbit_points()
    point_values = [
        tuple(
            sum(
                point[player]
                for player in range(N)
                if coalition >> player & 1
            )
            for coalition in range(GRAND + 1)
        )
        for point in points
    ]
    games = []
    game_set = set()
    best_margin = None
    edges = set()
    for batch in range(1, args.batches + 1):
        sample = rng.sample(range(len(points)), args.point_count)
        new_games = games_from_sample(point_values, sample)
        if args.permutation_close:
            new_games = {
                permute_game(game, permutation)
                for game in new_games
                for permutation in itertools.permutations(range(N))
            }
        for game in new_games:
            if game not in game_set:
                game_set.add(game)
                games.append(game)
        edges = set()
        add_induced_edges(games, edges)
        components = connected_components(len(games), edges)
        print(
            json.dumps(
                {
                    "batch": batch,
                    "phase": "solve",
                    "games": len(games),
                    "edges": len(edges),
                    "components": len(components),
                    "largest_component": len(components[0][0]) if components else 0,
                }
            ),
            flush=True,
        )
        best_margin = None
        for nodes, local_edges in components[:20]:
            margin = sparse_family_slack_float(
                [games[node] for node in nodes], local_edges, N, args.method
            )
            if margin is not None and (
                best_margin is None or margin < best_margin
            ):
                best_margin = margin
            if margin is not None and margin < -1e-8:
                break
        print(
            json.dumps(
                {
                    "batch": batch,
                    "games": len(games),
                    "edges": len(edges),
                    "components": len(components),
                    "largest_component": len(components[0][0]) if components else 0,
                    "margin": best_margin,
                }
            ),
            flush=True,
        )
        if best_margin is not None and best_margin < -1e-8:
            break
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best_margin is not None and best_margin < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "game_count": len(games),
        "edge_count": len(edges),
        "best_margin": best_margin,
        "games": [[str(value) for value in game] for game in games],
        "edges": [list(edge) for edge in sorted(edges)],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
