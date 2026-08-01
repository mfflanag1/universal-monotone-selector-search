#!/usr/bin/env python3
"""Greedily add permutation images using the sparse family LP."""

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
    parser.add_argument("--rounds", type=int, default=3)
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
    identity_image = tuple(sorted(base_games))
    images = {}
    for permutation in permutations(range(n)):
        image = tuple(
            sorted(
                permute_game(game, permutation)
                for game in base_games
            )
        )
        images.setdefault(image, permutation)
    remaining = {
        image: permutation
        for image, permutation in images.items()
        if image != identity_image
    }
    games = list(base_games)
    family = induced_family(
        games, n, "n5_sparse_permutation_greedy"
    )
    margin = sparse_family_slack_float(
        family.games, family.edges, n, args.method
    )
    history = []
    selected_permutations = [list(range(n))]
    for round_index in range(1, args.rounds + 1):
        best = None
        for image, permutation in remaining.items():
            candidate = induced_family(
                games + list(image),
                n,
                "n5_sparse_permutation_greedy",
            )
            candidate_margin = sparse_family_slack_float(
                candidate.games,
                candidate.edges,
                n,
                args.method,
            )
            if candidate_margin is None:
                raise RuntimeError("sparse greedy LP failed")
            if best is None or candidate_margin < best[0]:
                best = (
                    candidate_margin,
                    permutation,
                    image,
                    candidate,
                )
        if best is None:
            break
        margin, permutation, image, family = best
        games = family.games
        remaining.pop(image)
        selected_permutations.append(list(permutation))
        common_gap = common_core_gap_float(games, n)
        record = {
            "round": round_index,
            "permutation": list(permutation),
            "node_count": len(games),
            "edge_count": len(family.edges),
            "common_core_budget_gap_float": common_gap,
            "max_min_monotonicity_margin_float": margin,
        }
        history.append(record)
        print(json.dumps(record), flush=True)
        if margin < -1e-9:
            break

    common_gap = common_core_gap_float(games, n)
    result = {
        "status": "sparse_permutation_greedy_complete",
        "n": n,
        "source": str(args.input),
        "method": args.method,
        "distinct_images": len(images),
        "selected_permutations": selected_permutations,
        "rounds_completed": len(history),
        "node_count": len(games),
        "edge_count": len(family.edges),
        "common_core_budget_gap_float": common_gap,
        "max_min_monotonicity_margin_float": margin,
        "history": history,
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
                    "rounds_completed",
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
