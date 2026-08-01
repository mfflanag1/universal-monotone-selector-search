#!/usr/bin/env python3
"""Accumulate one-coordinate comparisons across random point-generated exact games."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from mixed_point_family_fast import random_points
from n5_bppv_pool_union_search import games_from_sample
from n5_point_deletion_family_search import connected_components
from n5_point_deletion_walk_search import add_induced_edges
from n5_sparse_family_lp import sparse_family_slack_float


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batches", type=int, default=100)
    parser.add_argument("--pool-count", type=int, default=24)
    parser.add_argument("--sample-count", type=int, default=9)
    parser.add_argument("--total", type=int, default=9)
    parser.add_argument("--bound", type=int, default=7)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--component-limit", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    n = 5
    grand = (1 << n) - 1
    points = random_points(rng, n, args.pool_count, args.total, args.bound)
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
    games: list[tuple[int, ...]] = []
    game_set: set[tuple[int, ...]] = set()
    edges: set[tuple[int, int, int]] = set()
    best = None
    for batch in range(1, args.batches + 1):
        sample = rng.sample(range(len(points)), args.sample_count)
        for game in games_from_sample(point_values, sample):
            if game not in game_set:
                game_set.add(game)
                games.append(game)
        edges.clear()
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
        batch_best = None
        for nodes, local_edges in components[: args.component_limit]:
            margin = sparse_family_slack_float(
                [games[node] for node in nodes], local_edges, n, "highs-ipm"
            )
            if margin is None:
                continue
            if batch_best is None or margin < batch_best:
                batch_best = margin
            if best is None or margin < best[0]:
                best = (
                    margin,
                    [games[node] for node in nodes],
                    list(local_edges),
                )
            if margin < -1e-8:
                break
        print(
            json.dumps(
                {
                    "batch": batch,
                    "games": len(games),
                    "edges": len(edges),
                    "components": len(components),
                    "largest_component": len(components[0][0]) if components else 0,
                    "margin": batch_best,
                }
            ),
            flush=True,
        )
        if batch_best is not None and batch_best < -1e-8:
            break
    found = best is not None and best[0] < -1e-8
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if found
            else "no_incompatible_exact_family_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "points": [list(point) for point in points],
        "game_count": len(games),
        "edge_count": len(edges),
        "best_margin": None if best is None else best[0],
        "games": [list(game) for game in best[1]] if found else None,
        "edges": [list(edge) for edge in best[2]] if found else None,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"status": payload["status"], "best_margin": payload["best_margin"]}))
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
