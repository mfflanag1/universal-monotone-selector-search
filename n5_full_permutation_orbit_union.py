#!/usr/bin/env python3
"""Solve the induced union of every player-permuted family image."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import permutations
from pathlib import Path

from family_search import (
    box_family_slack_float,
    common_core_gap_float,
    family_slack_float,
    serial_family,
)
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
    parser.add_argument("--skip-box", action="store_true")
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
    images = {
        tuple(
            sorted(
                permute_game(game, permutation)
                for game in base_games
            )
        )
        for permutation in permutations(range(n))
    }
    family = induced_family(
        [
            game
            for image in images
            for game in image
        ],
        n,
        "n5_full_permutation_orbit_union",
    )
    margin = sparse_family_slack_float(
        family.games, family.edges, n, args.method
    )
    if margin is None:
        raise RuntimeError("sparse permutation-orbit LP failed")
    box_margin = None
    if not args.skip_box:
        box_margin, _ = box_family_slack_float(
            family.games, family.edges, n
        )
    common_gap = common_core_gap_float(family.games, n)
    result = {
        "status": "full_permutation_orbit_union_complete",
        "n": n,
        "source": str(args.input),
        "distinct_images": len(images),
        "node_count": len(family.games),
        "edge_count": len(family.edges),
        "common_core_budget_gap_float": common_gap,
        "max_min_monotonicity_margin_float": margin,
        "box_max_min_monotonicity_margin_float": box_margin,
        "non_atomic_facet_tax_float": (
            box_margin - margin
            if box_margin is not None
            else None
        ),
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
                    "node_count",
                    "edge_count",
                    "common_core_budget_gap_float",
                    "max_min_monotonicity_margin_float",
                    "box_max_min_monotonicity_margin_float",
                    "non_atomic_facet_tax_float",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
