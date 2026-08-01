#!/usr/bin/env python3
"""Search asymmetric punctured cubes in the exact five-player cone."""

from __future__ import annotations

import argparse
import heapq
import json
import random
from pathlib import Path
from typing import Any

from n5_facet_search import exact_archive, search
from n5_topology_batch import punctured_cube_edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=100)
    parser.add_argument("--starts", type=int, default=20)
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--base-bump", type=float, default=0.02)
    parser.add_argument("--min-directions", type=int, default=3)
    parser.add_argument("--max-directions", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--retain", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pool = list(range(1, 32))
    summaries = []
    retained: list[tuple[float, int, dict[str, Any], list[int], list[float]]] = []
    for case in range(args.cases):
        count = rng.randint(args.min_directions, args.max_directions)
        directions = rng.sample(pool, count)
        direction_bumps = [
            args.base_bump * rng.randint(1, 8) / 2.0
            for _ in directions
        ]
        edges = punctured_cube_edges(directions)
        bump_by_coalition = dict(zip(directions, direction_bumps, strict=True))
        bumps = [bump_by_coalition[coalition] for _, _, coalition in edges]
        try:
            record = search(
                edges,
                bumps,
                args.starts,
                args.iterations,
                args.seed + case,
                quiet=True,
            )
        except RuntimeError:
            continue
        margin = float(record["best_margin_float"])
        summary = {
            "case": case,
            "directions": directions,
            "direction_bumps": direction_bumps,
            "node_count": 2**count - 1,
            "margin": margin,
        }
        summaries.append(summary)
        item = (-margin, case, record, directions, direction_bumps)
        if len(retained) < args.retain:
            heapq.heappush(retained, item)
        elif item > retained[0]:
            heapq.heapreplace(retained, item)
        print(json.dumps(summary), flush=True)
        if margin < -1e-8:
            break
    exact_cases = []
    for _, case, record, directions, direction_bumps in sorted(
        retained, key=lambda item: float(item[2]["best_margin_float"])
    ):
        archive = exact_archive(record, f"n5_random_punctured_cube_{case}")
        exact_cases.append(
            {
                "case": case,
                "directions": directions,
                "direction_bumps": direction_bumps,
                "archive": archive,
            }
        )
    payload = {
        "status": "random_exact_punctured_cube_complete",
        "configuration": vars(args) | {"output": str(args.output)},
        "summaries": summaries,
        "exact_cases": exact_cases,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
