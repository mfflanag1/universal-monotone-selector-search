#!/usr/bin/env python3
"""Search complete local exact-boundary stars around random exact games."""

from __future__ import annotations

import argparse
import json
import random
from fractions import Fraction
from pathlib import Path

from family_search import family_slack_float, game_from_points
from n5_exactified_obstruction_paths import exact_root
from n5_facet_search import load_facets
from n5_local_boundary_pool import coordinate_limit
from n5_point_deletion_family_search import compositions


F = Fraction
N = 5
GRAND = 31


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--point-count", type=int, default=8)
    parser.add_argument("--total", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pool = list(compositions(args.total))
    facets = load_facets()[0]
    rows = []
    best = None
    for trial in range(args.trials):
        central = game_from_points(rng.sample(pool, args.point_count))
        games = [central]
        edges = []
        for coalition in range(1, GRAND):
            for increase in (False, True):
                delta = coordinate_limit(
                    central, coalition, facets, increase
                )
                if delta <= 0:
                    continue
                neighbor = list(central)
                neighbor[coalition] += delta if increase else -delta
                games.append(tuple(neighbor))
                edge = (
                    (0, len(games) - 1, coalition)
                    if increase
                    else (len(games) - 1, 0, coalition)
                )
                edges.append(edge)
        root = exact_root(central, facets)
        if root is not None and root < central[GRAND]:
            neighbor = list(central)
            neighbor[GRAND] = root
            games.append(tuple(neighbor))
            edges.append((len(games) - 1, 0, GRAND))
        successor = list(central)
        successor[GRAND] += 1
        games.append(tuple(successor))
        edges.append((0, len(games) - 1, GRAND))
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
            "games": [[str(value) for value in game] for game in best[1]],
            "edges": [list(edge) for edge in best[2]],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
