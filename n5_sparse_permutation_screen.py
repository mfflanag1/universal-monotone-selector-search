#!/usr/bin/env python3
"""Screen permutation pairs with the memory-safe sparse family LP."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import permutations
from pathlib import Path

from family_search import common_core_gap_float, serial_family
from n5_intrinsic_permutation_union import (
    induced_family,
    permute_game,
)
from n5_sparse_family_lp import sparse_family_slack_float


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--method",
        choices=("highs", "highs-ds", "highs-ipm"),
        default="highs-ipm",
    )
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    n = int(payload["n"])
    base_games = [
        tuple(F(value) for value in game)
        for game in payload["family"]["games"]
    ]
    identity = tuple(range(n))
    images = {}
    for permutation in permutations(range(n)):
        if permutation == identity:
            continue
        image = tuple(
            sorted(
                permute_game(game, permutation)
                for game in base_games
            )
        )
        images.setdefault(image, permutation)

    rows = []
    best = None
    for index, (image, permutation) in enumerate(
        images.items(), start=1
    ):
        family = induced_family(
            base_games + list(image),
            n,
            "n5_sparse_permutation_pair",
        )
        margin = sparse_family_slack_float(
            family.games, family.edges, n, args.method
        )
        if margin is None:
            raise RuntimeError("sparse permutation LP failed")
        row = {
            "permutation": list(permutation),
            "node_count": len(family.games),
            "edge_count": len(family.edges),
            "max_min_monotonicity_margin_float": margin,
        }
        rows.append(row)
        if best is None or margin < best[0]:
            best = (margin, permutation, family)
            print(
                json.dumps(
                    {
                        "screened": index,
                        "distinct_images": len(images),
                        **row,
                    }
                ),
                flush=True,
            )

    if best is None:
        raise RuntimeError("no distinct permutation image")
    margin, permutation, family = best
    common_gap = common_core_gap_float(family.games, n)
    result = {
        "status": "sparse_permutation_pair_screen_complete",
        "n": n,
        "source": str(args.input),
        "method": args.method,
        "distinct_images": len(images),
        "best_permutation": list(permutation),
        "node_count": len(family.games),
        "edge_count": len(family.edges),
        "common_core_budget_gap_float": common_gap,
        "max_min_monotonicity_margin_float": margin,
        "rows": rows,
        "family": serial_family(family, margin, common_gap),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "distinct_images",
                    "best_permutation",
                    "node_count",
                    "edge_count",
                    "common_core_budget_gap_float",
                    "max_min_monotonicity_margin_float",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
