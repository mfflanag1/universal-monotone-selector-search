#!/usr/bin/env python3
"""Search aggregate caps with complete local exact stars at six players."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np

from family_search import family_slack_float
from n6_local_closure_search import moves
from n6_local_star_search import N, point_game
from n6_mixed_terminal_search import load_facets


def random_points(rng: random.Random, count: int, bound: int):
    translation = 5 * bound
    points = []
    while len(points) < count:
        prefix = [rng.randint(-bound, bound) for _ in range(N - 1)]
        point = tuple(
            value + translation for value in (*prefix, -sum(prefix))
        )
        if point not in points:
            points.append(point)
    return points


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--point-count", type=int, default=10)
    parser.add_argument("--coordinate-bound", type=int, default=5)
    parser.add_argument("--grand-bump", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    facets = np.asarray(load_facets(), dtype=float)
    rows = []
    best = None
    for trial in range(args.trials):
        lower = point_game(
            random_points(
                rng, args.point_count, args.coordinate_bound
            )
        )
        upper = lower.copy()
        upper[63] += args.grand_bump
        if float((facets @ upper).min()) < -1e-7:
            continue
        games = [lower, upper]
        edges = [(0, 1, 63)]
        for neighbor, coalition, increase in moves(upper, facets):
            target = len(games)
            games.append(neighbor)
            edges.append(
                (1, target, coalition)
                if increase
                else (target, 1, coalition)
            )
        margin, _ = family_slack_float(games, edges, N)
        row = {
            "trial": trial,
            "nodes": len(games),
            "edges": len(edges),
            "margin": margin,
        }
        rows.append(row)
        if best is None or margin < best[0]:
            best = (margin, games, edges)
            print(json.dumps(row), flush=True)
        if margin < -1e-8:
            break
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "rows": rows,
        "best": None
        if best is None
        else {
            "margin": best[0],
            "games": [[float(value) for value in game] for game in best[1]],
            "edges": [list(edge) for edge in best[2]],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
