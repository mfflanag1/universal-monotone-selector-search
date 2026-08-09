#!/usr/bin/env python3
"""Search degree-three K3,3 protected block flows on five players."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any

from n5_mixed_terminal_search import complete_bipartite_topology, optimize_topology


N = 5
GRAND = (1 << N) - 1


def draw_candidate(rng: random.Random) -> tuple[int, ...]:
    while True:
        masks = tuple(rng.sample(range(1, GRAND), 9))
        topology = complete_bipartite_topology(masks, (1,) * 9, 3, N)
        level_counts = [len(set(vector)) for vector in topology.divergence]
        if min(level_counts) < 3:
            continue
        if sum(level_counts) < 21:
            continue
        if masks[0] & masks[1] & masks[2] & masks[3] & masks[4] & masks[5] & masks[6] & masks[7] & masks[8]:
            continue
        return masks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260816)
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
        topology = complete_bipartite_topology(masks, (1,) * 9, 3, N)
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
            "terminal_level_counts": [
                len(set(vector)) for vector in topology.divergence
            ],
            "divergence": [list(vector) for vector in topology.divergence],
            **optimum,
        }
        if float(row["objective"]) > 1e-8:
            positive += 1
        best.append(row)
        best.sort(key=lambda item: float(item["objective"]), reverse=True)
        del best[10:]
        if len(seen) % 10 == 0 or positive:
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
        "topology": "complete_bipartite_k33_high_terminal_diversity",
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
