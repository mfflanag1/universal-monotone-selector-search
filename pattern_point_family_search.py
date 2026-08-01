#!/usr/bin/env python3
"""Search subset-generated exact families from a repeated-coordinate orbit."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from fractions import Fraction
from pathlib import Path

import n5_point_deletion_family_search as deletion
from n5_sparse_family_lp import sparse_family_slack_float


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", required=True)
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--point-count", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pattern = tuple(F(value) for value in args.pattern.split(","))
    n = len(pattern)
    deletion.N = n
    pool = sorted(
        {
            tuple(pattern[index] for index in permutation)
            for permutation in itertools.permutations(range(n))
        }
    )
    if args.point_count > len(pool):
        raise ValueError("point count exceeds orbit size")
    rng = random.Random(args.seed)
    rows = []
    best = None
    for trial in range(args.trials):
        points = rng.sample(pool, args.point_count)
        games, edges = deletion.build(points)
        components = deletion.connected_components(len(games), edges)
        trial_best = None
        for nodes, local_edges in components[:10]:
            margin = sparse_family_slack_float(
                [games[node] for node in nodes], local_edges, n, "highs-ipm"
            )
            if margin is None:
                continue
            if trial_best is None or margin < trial_best[0]:
                trial_best = (margin, nodes, local_edges)
            if best is None or margin < best[0]:
                best = (
                    margin,
                    points,
                    [games[node] for node in nodes],
                    local_edges,
                )
            if margin < -1e-8:
                break
        row = {
            "trial": trial,
            "games": len(games),
            "edges": len(edges),
            "largest_component": len(components[0][0]) if components else 0,
            "margin": None if trial_best is None else trial_best[0],
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
        if trial_best is not None and trial_best[0] < -1e-8:
            break
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "pattern": [str(value) for value in pattern],
        "orbit_size": len(pool),
        "rows": rows,
        "best": None
        if best is None
        else {
            "margin": best[0],
            "points": [[str(value) for value in point] for point in best[1]],
            "games": [[str(value) for value in game] for game in best[2]],
            "edges": [list(edge) for edge in best[3]],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
