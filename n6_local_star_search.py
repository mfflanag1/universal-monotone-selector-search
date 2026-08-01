#!/usr/bin/env python3
"""Search complete proper-coordinate boundary stars in the n=6 exact cone."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from pathlib import Path

import numpy as np

from family_search import family_slack_float
from n6_mixed_terminal_search import load_facets


N = 6
GRAND = 63


def compositions(total: int):
    for cuts in itertools.combinations(range(total + N - 1), N - 1):
        values = []
        previous = -1
        for cut in (*cuts, total + N - 1):
            values.append(cut - previous - 1)
            previous = cut
        yield tuple(values)


def point_game(points: list[tuple[int, ...]]) -> np.ndarray:
    return np.asarray(
        [
            min(
                sum(point[player] for player in range(N) if mask >> player & 1)
                for point in points
            )
            for mask in range(64)
        ],
        dtype=float,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--point-count", type=int, default=10)
    parser.add_argument("--total", type=int, default=10)
    parser.add_argument("--coordinate-bound", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pool = list(compositions(args.total)) if not args.coordinate_bound else []
    facets = np.asarray(load_facets(), dtype=float)
    rows = []
    best = None
    for trial in range(args.trials):
        if args.coordinate_bound:
            translation = 5 * args.coordinate_bound
            points = []
            while len(points) < args.point_count:
                prefix = [
                    rng.randint(-args.coordinate_bound, args.coordinate_bound)
                    for _ in range(N - 1)
                ]
                point = tuple(
                    value + translation
                    for value in (*prefix, -sum(prefix))
                )
                if point not in points:
                    points.append(point)
            central = point_game(points)
        else:
            central = point_game(rng.sample(pool, args.point_count))
        slack = facets @ central
        if float(slack.min()) < -1e-7:
            raise RuntimeError("point game is outside the exact cone")
        games = [central.copy()]
        edges = []
        for coalition in range(1, GRAND):
            increase_bounds = []
            negative = facets[:, coalition] < 0
            if negative.any():
                increase_bounds.append(
                    float(np.min(slack[negative] / -facets[negative, coalition]))
                )
            increase_bounds.extend(
                central[upper] - central[coalition]
                for upper in range(1, 64)
                if upper != coalition and upper & coalition == coalition
            )
            increase = min(increase_bounds)
            if increase > 1e-8:
                neighbor = central.copy()
                neighbor[coalition] += increase
                games.append(neighbor)
                edges.append((0, len(games) - 1, coalition))

            decrease_bounds = []
            positive = facets[:, coalition] > 0
            if positive.any():
                decrease_bounds.append(
                    float(np.min(slack[positive] / facets[positive, coalition]))
                )
            decrease_bounds.extend(
                central[coalition] - central[lower]
                for lower in range(0, GRAND)
                if lower != coalition and lower & coalition == lower
            )
            decrease = min(decrease_bounds)
            if decrease > 1e-8:
                neighbor = central.copy()
                neighbor[coalition] -= decrease
                games.append(neighbor)
                edges.append((len(games) - 1, 0, coalition))
        successor = central.copy()
        successor[GRAND] += 1.0
        if float((facets @ successor).min()) >= -1e-7:
            games.append(successor)
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
            "games": [[float(value) for value in game] for game in best[1]],
            "edges": [list(edge) for edge in best[2]],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
