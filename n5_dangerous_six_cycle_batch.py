#!/usr/bin/env python3
"""Search six-cycles whose six terminal divergences all have sharp premium 1/2."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any

from n5_mixed_terminal_search import (
    bipartite_six_cycle_topology,
    optimize_topology,
)


N = 5
GRAND = (1 << N) - 1
MAXIMUM_PREMIUM_TYPES = {(1, 2, 2), (1, 3, 1), (2, 2, 1)}


def level_multiplicities(left: int, right: int, negative: bool) -> tuple[int, int, int]:
    intersection = (left & right).bit_count()
    symmetric_difference = (left ^ right).bit_count()
    neither = N - (left | right).bit_count()
    return (
        (intersection, symmetric_difference, neither)
        if negative
        else (neither, symmetric_difference, intersection)
    )


def all_terminals_sharp(masks: tuple[int, ...]) -> bool:
    # Source pairs at nodes 0, 1, 2; sink pairs at nodes 3, 4, 5.
    signed_pairs = (
        (masks[0], masks[5], False),
        (masks[1], masks[2], False),
        (masks[3], masks[4], False),
        (masks[0], masks[1], True),
        (masks[2], masks[3], True),
        (masks[4], masks[5], True),
    )
    return all(
        level_multiplicities(left, right, negative) in MAXIMUM_PREMIUM_TYPES
        for left, right, negative in signed_pairs
    )


def draw_candidate(rng: random.Random) -> tuple[int, ...]:
    while True:
        masks = tuple(rng.sample(range(1, GRAND), 6))
        if not all_terminals_sharp(masks):
            continue
        if masks[0] & masks[1] & masks[2] & masks[3] & masks[4] & masks[5]:
            continue
        if masks[0] | masks[1] | masks[2] | masks[3] | masks[4] | masks[5] != GRAND:
            continue
        return masks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument("--iterations", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--variable-grand", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    seen: set[tuple[int, ...]] = set()
    best: list[dict[str, Any]] = []
    positive = 0
    while len(seen) < args.samples:
        masks = draw_candidate(rng)
        if masks in seen:
            continue
        seen.add(masks)
        topology = bipartite_six_cycle_topology(masks, N)
        optimum = optimize_topology(
            topology,
            N,
            rng,
            args.starts,
            args.iterations,
            args.variable_grand,
        )
        if optimum is None:
            continue
        row = {
            "masks": list(masks),
            "divergence": [list(vector) for vector in topology.divergence],
            **optimum,
        }
        if float(row["objective"]) > 1e-8:
            positive += 1
        best.append(row)
        best.sort(key=lambda item: float(item["objective"]), reverse=True)
        del best[10:]
        if len(seen) % 25 == 0 or positive:
            print(
                json.dumps(
                    {
                        "tested": len(seen),
                        "positive": positive,
                        "best_objective": best[0]["objective"],
                    }
                ),
                flush=True,
            )
        if positive:
            break
    payload = {
        "status": (
            "positive_terminal_obstruction_found"
            if positive
            else "no_positive_terminal_obstruction_found"
        ),
        "n": N,
        "topology": "alternating_six_cycle_all_terminals_sharp",
        "tested": len(seen),
        "positive": positive,
        "variable_grand": args.variable_grand,
        "configuration": vars(args) | {"output": str(args.output)},
        "best": best,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 1 if positive else 0


if __name__ == "__main__":
    raise SystemExit(main())
