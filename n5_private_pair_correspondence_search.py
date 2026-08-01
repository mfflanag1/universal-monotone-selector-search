#!/usr/bin/env python3
"""Search private-facet exact pairs for failure of core-correspondence monotonicity."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np

from n5_core_correspondence_edge_audit import extendable
from private_facet_search import game_from_points, generate_points, private_minimizer


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pairs = 0
    witness = None
    for trial in range(1, args.trials + 1):
        points = generate_points(5, rng.randint(3, 12), rng.randint(5, 50), 20, 0.7, rng)
        upper = game_from_points(points)
        coalitions = list(range(1, 31))
        rng.shuffle(coalitions)
        for coalition in coalitions:
            private = private_minimizer(5, upper, coalition)
            if private is None:
                continue
            lower_value, point = private
            if upper[coalition] - lower_value <= 1e-7:
                continue
            pairs += 1
            lower = list(upper)
            lower[coalition] = lower_value
            if not extendable(upper, np.asarray(point), coalition):
                witness = {
                    "trial": trial,
                    "coalition": coalition,
                    "generators": points,
                    "private_point": point,
                    "lower_game": lower,
                    "upper_game": upper,
                }
                break
        if witness is not None:
            break
        if trial % 100 == 0:
            print(json.dumps({"trials": trial, "pairs": pairs}), flush=True)
    payload = {
        "status": "nonextendable_private_point_found" if witness else "none_found",
        "trials": trial,
        "pairs": pairs,
        "witness": witness,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "trials", "pairs")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
