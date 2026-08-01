#!/usr/bin/env python3
"""Fast exact-family search over asymmetric equal-total witness points."""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_point_deletion_family_search import connected_components
from n5_sparse_family_lp import sparse_family_slack_float


def build(points: list[tuple[int, ...]], induced_edges: bool = False):
    n = len(points[0])
    grand = (1 << n) - 1
    point_values = np.asarray(
        [
            [
                sum(point[player] for player in range(n) if coalition >> player & 1)
                for coalition in range(grand + 1)
            ]
            for point in points
        ],
        dtype=np.int16,
    )
    state_count = 1 << len(points)
    state_games = np.full((state_count, grand + 1), 32767, dtype=np.int16)
    for state in range(1, state_count):
        bit = state & -state
        point = bit.bit_length() - 1
        state_games[state] = np.minimum(state_games[state ^ bit], point_values[point])
    games, inverse = np.unique(state_games[1:], axis=0, return_inverse=True)
    state_game = np.empty(state_count, dtype=np.int32)
    state_game[0] = -1
    state_game[1:] = inverse
    edges: set[tuple[int, int, int]] = set()
    for state in range(1, state_count):
        if state & (state - 1) == 0:
            continue
        lower = int(state_game[state])
        remaining = state
        while remaining:
            bit = remaining & -remaining
            upper = int(state_game[state ^ bit])
            if lower != upper:
                differences = np.flatnonzero(games[lower] != games[upper])
                if len(differences) == 1:
                    coalition = int(differences[0])
                    if 0 < coalition < grand and games[lower, coalition] < games[upper, coalition]:
                        edges.add((lower, upper, coalition))
            remaining ^= bit
    if induced_edges:
        for coalition in range(1, grand):
            groups: dict[bytes, list[int]] = {}
            keep = [mask for mask in range(grand + 1) if mask != coalition]
            for game_index, game in enumerate(games):
                groups.setdefault(game[keep].tobytes(), []).append(game_index)
            for indices in groups.values():
                ordered = sorted(indices, key=lambda index: int(games[index, coalition]))
                for position, lower in enumerate(ordered):
                    for upper in ordered[position + 1 :]:
                        if games[lower, coalition] < games[upper, coalition]:
                            edges.add((lower, upper, coalition))
    return games, edges


def random_points(
    rng: random.Random, n: int, count: int, total: int, bound: int
) -> list[tuple[int, ...]]:
    points: set[tuple[int, ...]] = set()
    while len(points) < count:
        prefix = tuple(rng.randint(-bound, bound) for _ in range(n - 1))
        last = total - sum(prefix)
        if -2 * bound <= last <= 2 * bound:
            points.add(prefix + (last,))
    return list(points)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--point-count", type=int, default=12)
    parser.add_argument("--total", type=int, default=10)
    parser.add_argument("--bound", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--induced-edges", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    rows = []
    best = None
    for trial in range(args.trials):
        points = random_points(
            rng, args.n, args.point_count, args.total, args.bound
        )
        games_array, edges = build(points, args.induced_edges)
        components = connected_components(len(games_array), edges)
        trial_best = None
        for nodes, local_edges in components[:20]:
            games = [tuple(int(value) for value in games_array[node]) for node in nodes]
            margin = sparse_family_slack_float(games, local_edges, args.n, "highs-ipm")
            if margin is None:
                continue
            if trial_best is None or margin < trial_best:
                trial_best = margin
            if best is None or margin < best[0]:
                best = (margin, points, games, local_edges)
            if margin < -1e-8:
                break
        row = {
            "trial": trial,
            "games": len(games_array),
            "edges": len(edges),
            "components": len(components),
            "largest_component": len(components[0][0]) if components else 0,
            "margin": trial_best,
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
        if trial_best is not None and trial_best < -1e-8:
            break
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "rows": rows,
        "best": None
        if best is None
        else {
            "margin": best[0],
            "points": [list(point) for point in best[1]],
            "games": [list(game) for game in best[2]] if best[0] < -1e-8 else None,
            "edges": [list(edge) for edge in best[3]] if best[0] < -1e-8 else None,
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"status": payload["status"], "best_margin": None if best is None else best[0]}))
    return 1 if best is not None and best[0] < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
