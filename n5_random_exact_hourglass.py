#!/usr/bin/env python3
"""Adversarially optimize random exact five-player hourglass families."""

from __future__ import annotations

import argparse
import heapq
import json
import random
from pathlib import Path
from typing import Any

from n5_facet_search import exact_archive, search


def coalition_pool(full: bool) -> list[int]:
    if full:
        return list(range(1, 32))
    return [mask for mask in range(1, 31) if 2 <= mask.bit_count() <= 4]


def hourglass(
    incoming: list[int], outgoing: list[int]
) -> list[tuple[int, int, int]]:
    edges = [
        (index + 1, 0, coalition)
        for index, coalition in enumerate(incoming)
    ]
    offset = 1 + len(incoming)
    edges.extend(
        (0, offset + index, coalition)
        for index, coalition in enumerate(outgoing)
    )
    return edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=100)
    parser.add_argument("--starts", type=int, default=20)
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--bump", type=float, default=0.05)
    parser.add_argument("--asymmetric", action="store_true")
    parser.add_argument("--full-pool", action="store_true")
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--retain", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pool = coalition_pool(args.full_pool)
    retained: list[tuple[float, int, dict[str, Any], list[int], list[int]]] = []
    summaries = []
    for case in range(args.cases):
        incoming = rng.sample(pool, rng.randint(2, 6))
        outgoing = rng.sample(pool, rng.randint(2, 6))
        edges = hourglass(incoming, outgoing)
        bumps = [
            args.bump * (rng.randint(1, 12) / 4.0)
            if args.asymmetric
            else args.bump
            for _ in edges
        ]
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
            "incoming": incoming,
            "outgoing": outgoing,
            "bumps": bumps,
            "margin": margin,
        }
        summaries.append(summary)
        item = (-margin, case, record, incoming, outgoing)
        if len(retained) < args.retain:
            heapq.heappush(retained, item)
        elif item > retained[0]:
            heapq.heapreplace(retained, item)
        print(json.dumps(summary), flush=True)
        if margin < -1e-8:
            break
    exact_cases = []
    for _, case, record, incoming, outgoing in sorted(
        retained, key=lambda item: float(item[2]["best_margin_float"])
    ):
        archive = exact_archive(record, f"n5_random_hourglass_{case}")
        exact_cases.append(
            {
                "case": case,
                "incoming": incoming,
                "outgoing": outgoing,
                "archive": archive,
            }
        )
    payload = {
        "status": "random_exact_hourglass_complete",
        "configuration": {
            "cases": args.cases,
            "starts": args.starts,
            "iterations": args.iterations,
            "bump": args.bump,
            "asymmetric": args.asymmetric,
            "full_pool": args.full_pool,
            "seed": args.seed,
        },
        "summaries": summaries,
        "exact_cases": exact_cases,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
