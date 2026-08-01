#!/usr/bin/env python3
"""Search six-player private-facet pairs for a nonextendable lower core point."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from private_facet_search import game_from_points, generate_points, private_minimizer


def extendable(game: list[float], point: list[float], protected: int) -> bool:
    n = 6
    grand = (1 << n) - 1
    rows = []
    rhs = []
    for coalition in range(1, grand):
        rows.append([-float(bool(coalition >> player & 1)) for player in range(n)])
        rhs.append(-float(game[coalition]))
    for player in range(n):
        if protected >> player & 1:
            row = [0.0] * n
            row[player] = -1.0
            rows.append(row)
            rhs.append(-float(point[player]))
    result = linprog(
        np.zeros(n),
        A_ub=np.asarray(rows),
        b_ub=np.asarray(rhs),
        A_eq=np.ones((1, n)),
        b_eq=[float(game[grand])],
        bounds=[(None, None)] * n,
        method="highs",
    )
    return result.success


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pairs = 0
    witness = None
    five_player_coalitions = [63 ^ (1 << player) for player in range(6)]
    for trial in range(1, args.trials + 1):
        points = generate_points(6, rng.randint(3, 14), rng.randint(6, 60), 25, 0.8, rng)
        upper = game_from_points(points)
        rng.shuffle(five_player_coalitions)
        for coalition in five_player_coalitions:
            private = private_minimizer(6, upper, coalition)
            if private is None:
                continue
            lower_value, point = private
            if upper[coalition] - lower_value <= 1e-7:
                continue
            pairs += 1
            lower = list(upper)
            lower[coalition] = lower_value
            if not extendable(upper, point, coalition):
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
