#!/usr/bin/env python3
"""Adversarial exact-cone search on six-player aggregate-capped bowties."""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path


sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n6_mixed_terminal_search import load_facets

from n6_exact_mixed_cube_search import GRAND, optimize


def capped_bowtie(directions):
    edges = [(0, 1, GRAND)]
    for arm, (first, second) in enumerate(directions):
        left = 2 + 3 * arm
        right = left + 1
        top = left + 2
        edges.extend(
            (
                (1, left, first),
                (1, right, second),
                (left, top, second),
                (right, top, first),
            )
        )
    return edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=50)
    parser.add_argument("--arms", type=int, default=2)
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--bump", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facets = load_facets()
    rng = random.Random(args.seed)
    pool = [mask for mask in range(1, GRAND) if 2 <= mask.bit_count() <= 5]
    rows = []
    best = None
    for case in range(args.cases):
        directions = tuple(tuple(rng.sample(pool, 2)) for _ in range(args.arms))
        edges = capped_bowtie(directions)
        bump_by_coalition = {
            coalition: args.bump * rng.randint(1, 5)
            for coalition in {changed for _lower, _upper, changed in edges}
        }
        bumps = [bump_by_coalition[changed] for _lower, _upper, changed in edges]
        margin, games = optimize(
            edges,
            bumps,
            facets,
            args.starts,
            args.iterations,
            args.seed + case,
        )
        if games is None:
            continue
        row = {
            "case": case,
            "directions": [list(pair) for pair in directions],
            "edges": [list(edge) for edge in edges],
            "bumps": bumps,
            "margin": margin,
        }
        rows.append(row)
        if best is None or margin < best["margin"]:
            best = row | {"games": games}
            print(json.dumps({"best": row}), flush=True)
        if margin < -1e-8:
            break
    payload = {
        "status": (
            "negative_margin_found"
            if best is not None and best["margin"] < -1e-8
            else "none_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "best": best,
        "rows": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
