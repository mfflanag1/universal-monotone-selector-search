#!/usr/bin/env python3
"""Test whether protected exact relaxations change mixed upper expectations."""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from nonlarge_complex_search import biswas_nonlarge_game


N = 5
GRAND = (1 << N) - 1


def core_optimum(game: np.ndarray, objective: np.ndarray, maximize: bool) -> float:
    rows = []
    rhs = []
    for coalition in range(1, GRAND):
        row = np.zeros(N)
        for player in range(N):
            if coalition >> player & 1:
                row[player] = -1.0
        rows.append(row)
        rhs.append(-game[coalition])
    result = linprog(
        -objective if maximize else objective,
        A_ub=np.asarray(rows),
        b_ub=np.asarray(rhs),
        A_eq=np.ones((1, N)),
        b_eq=np.asarray([game[GRAND]]),
        bounds=[(None, None)] * N,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return float((-1.0 if maximize else 1.0) * result.fun)


def relaxed_core_optimum(
    game: np.ndarray, support: int, objective: np.ndarray, maximize: bool
) -> float:
    rows = []
    rhs = []
    for coalition in range(1, GRAND):
        if coalition & support == support:
            continue
        row = np.zeros(N)
        for player in range(N):
            if coalition >> player & 1:
                row[player] = -1.0
        rows.append(row)
        rhs.append(-game[coalition])
    result = linprog(
        -objective if maximize else objective,
        A_ub=np.asarray(rows),
        b_ub=np.asarray(rhs),
        A_eq=np.ones((1, N)),
        b_eq=np.asarray([game[GRAND]]),
        bounds=[(None, None)] * N,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return float((-1.0 if maximize else 1.0) * result.fun)


def exactification(game: np.ndarray, support: int) -> list[float]:
    values = [0.0] * (GRAND + 1)
    values[GRAND] = float(game[GRAND])
    for coalition in range(1, GRAND):
        indicator = np.asarray(
            [float(coalition >> player & 1) for player in range(N)]
        )
        values[coalition] = relaxed_core_optimum(
            game, support, indicator, False
        )
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weight-bound", type=int, default=3)
    parser.add_argument("--random-cases", type=int, default=100)
    parser.add_argument("--random-tests-per-case", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260801)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    game = np.asarray(
        [float(value) for value in biswas_nonlarge_game(Fraction(2))]
    )
    rows = []
    for support in range(1, GRAND):
        players = [player for player in range(N) if support >> player & 1]
        if len(players) < 2:
            continue
        for positive_weights in itertools.product(
            range(1, args.weight_bound + 1), repeat=len(players)
        ):
            if len(set(positive_weights)) == 1:
                continue
            objective = np.zeros(N)
            for player, weight in zip(players, positive_weights, strict=True):
                objective[player] = weight
            original = core_optimum(game, objective, True)
            relaxed = relaxed_core_optimum(
                game, support, objective, True
            )
            gap = relaxed - original
            rows.append(
                {
                    "source": "biswas",
                    "support": support,
                    "weights": [int(value) for value in objective],
                    "original_upper": original,
                    "relaxed_upper": relaxed,
                    "gap": gap,
                }
            )
    rng = random.Random(args.seed)
    for case in range(args.random_cases):
        points = []
        for _ in range(rng.randint(2, 8)):
            cuts = sorted(rng.sample(range(1, 30), N - 1))
            points.append(
                np.diff(np.asarray((0, *cuts, 30), dtype=float))
            )
        random_game = np.asarray(
            [
                min(
                    sum(
                        point[player]
                        for player in range(N)
                        if coalition >> player & 1
                    )
                    for point in points
                )
                for coalition in range(GRAND + 1)
            ]
        )
        for _ in range(args.random_tests_per_case):
            support_players = rng.sample(range(N), rng.randint(2, N - 1))
            support = sum(1 << player for player in support_players)
            positive_weights = [
                rng.randint(1, args.weight_bound) for _ in support_players
            ]
            if len(set(positive_weights)) == 1:
                continue
            objective = np.zeros(N)
            for player, weight in zip(
                support_players, positive_weights, strict=True
            ):
                objective[player] = weight
            original = core_optimum(random_game, objective, True)
            relaxed = relaxed_core_optimum(
                random_game, support, objective, True
            )
            rows.append(
                {
                    "source": f"random_{case}",
                    "support": support,
                    "weights": [int(value) for value in objective],
                    "original_upper": original,
                    "relaxed_upper": relaxed,
                    "gap": relaxed - original,
                    "source_game": [float(value) for value in random_game],
                }
            )
    rows.sort(key=lambda row: float(row["gap"]), reverse=True)
    best = rows[0]
    best["relaxed_exactification"] = exactification(
        np.asarray(best.get("source_game", game)), int(best["support"])
    )
    payload = {
        "status": (
            "mixed_upper_drop_available"
            if float(best["gap"]) > 1e-9
            else "mixed_upper_preserved"
        ),
        "tested": len(rows),
        "best": best,
        "top": rows[:20],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key != "top"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
