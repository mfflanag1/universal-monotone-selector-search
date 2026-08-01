#!/usr/bin/env python3
"""Search the complete proper-coalition star over the exact five-player cone."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from family_search import common_core_gap_float
from n5_facet_search import search


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bump", type=float, default=0.001)
    parser.add_argument("--starts", type=int, default=16)
    parser.add_argument("--iterations", type=int, default=24)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--include-grand", action="store_true")
    args = parser.parse_args()

    upper = 32 if args.include_grand else 31
    edges = [(0, coalition, coalition) for coalition in range(1, upper)]
    record = search(
        edges,
        [args.bump] * len(edges),
        args.starts,
        args.iterations,
        args.seed,
    )
    result = {
        "status": (
            "n5_complete_star_float"
            if args.include_grand
            else "n5_complete_proper_star_float"
        ),
        "common_core_gap_float": (
            None
            if args.include_grand
            else common_core_gap_float(record["best_games_float"], 5)
        ),
        "record": record,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "margin": record["best_margin_float"],
                "common_core_gap": result["common_core_gap_float"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
