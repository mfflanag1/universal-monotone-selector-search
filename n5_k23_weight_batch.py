#!/usr/bin/env python3
"""Search scalar block weights on a fixed five-player K2,3 topology."""

from __future__ import annotations

import argparse
import json
import random
from functools import reduce
from math import gcd
from pathlib import Path
from typing import Any

from n5_mixed_terminal_search import (
    Topology,
    complete_bipartite_rectangular_topology,
    optimize_topology,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--masks", type=int, nargs=6, required=True)
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--max-weight", type=int, default=6)
    parser.add_argument("--starts", type=int, default=8)
    parser.add_argument("--iterations", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260827)
    parser.add_argument("--variable-grand", action="store_true")
    parser.add_argument("--player-specific", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    masks = tuple(args.masks)
    seen: set[tuple[int, ...]] = set()
    best: list[dict[str, Any]] = []
    positive = 0
    while len(seen) < args.samples:
        if args.player_specific:
            raw = tuple(
                rng.randint(1, args.max_weight)
                for mask in masks
                for player in range(5)
                if mask >> player & 1
            )
        else:
            raw = tuple(rng.randint(1, args.max_weight) for _ in range(6))
        common = reduce(gcd, raw)
        normalized = tuple(weight // common for weight in raw)
        if normalized in seen:
            continue
        seen.add(normalized)
        if args.player_specific:
            arcs = tuple(
                (source, 2 + sink, masks[source * 3 + sink])
                for source in range(2)
                for sink in range(3)
            )
            divergence = [[0] * 5 for _ in range(5)]
            cursor = 0
            for lower, upper, mask in arcs:
                for player in range(5):
                    if not (mask >> player & 1):
                        continue
                    weight = normalized[cursor]
                    cursor += 1
                    divergence[lower][player] += weight
                    divergence[upper][player] -= weight
            topology = Topology(
                5,
                arcs,
                tuple(tuple(vector) for vector in divergence),
            )
        else:
            topology = complete_bipartite_rectangular_topology(
                masks, normalized, 2, 3, 5
            )
        optimum = optimize_topology(
            topology,
            5,
            rng,
            args.starts,
            args.iterations,
            args.variable_grand,
        )
        if optimum is None:
            continue
        row = {
            "weights": list(normalized),
            "player_specific": args.player_specific,
            "divergence": [list(vector) for vector in topology.divergence],
            **optimum,
        }
        if float(row["objective"]) > 1e-8:
            positive += 1
        best.append(row)
        best.sort(key=lambda item: float(item["objective"]), reverse=True)
        del best[20:]
        if len(seen) % 20 == 0 or positive:
            print(
                json.dumps(
                    {
                        "tested": len(seen),
                        "positive": positive,
                        "best_objective": best[0]["objective"],
                        "best_weights": best[0]["weights"],
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
        "masks": list(masks),
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
