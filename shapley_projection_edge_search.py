#!/usr/bin/env python3
"""Test Euclidean projection of the Shapley value onto exact-game cores."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_facet_search import load_facets

from alexia_exact_edge_search import maximum_bump
from family_search import game_from_points, random_composition


F = Fraction
N = 5
GRAND = (1 << N) - 1


def shapley(game: tuple[F, ...]) -> np.ndarray:
    result = np.zeros(N)
    denominator = math.factorial(N)
    for player in range(N):
        for coalition in range(GRAND + 1):
            if coalition >> player & 1:
                continue
            size = coalition.bit_count()
            weight = (
                math.factorial(size)
                * math.factorial(N - size - 1)
                / denominator
            )
            result[player] += weight * float(
                game[coalition | (1 << player)] - game[coalition]
            )
    return result


def projected_shapley(game: tuple[F, ...]) -> np.ndarray:
    reference = shapley(game)
    inequalities = []
    for coalition in range(1, GRAND):
        row = np.asarray(
            [float(coalition >> player & 1) for player in range(N)]
        )
        worth = float(game[coalition])
        inequalities.append(
            {
                "type": "ineq",
                "fun": lambda x, row=row, worth=worth: float(row @ x - worth),
                "jac": lambda _x, row=row: row,
            }
        )
    constraints = inequalities + [
        {
            "type": "eq",
            "fun": lambda x: float(np.sum(x) - float(game[GRAND])),
            "jac": lambda _x: np.ones(N),
        }
    ]
    start = np.asarray(
        [float(value) for value in next(iter(core_witnesses(game)))],
        dtype=float,
    )
    solved = minimize(
        lambda x: 0.5 * float(np.sum((x - reference) ** 2)),
        start,
        jac=lambda x: x - reference,
        constraints=constraints,
        method="SLSQP",
        options={"ftol": 1e-12, "maxiter": 1000},
    )
    if not solved.success or min(constraint["fun"](solved.x) for constraint in inequalities) < -1e-7:
        raise RuntimeError(solved.message)
    return solved.x


def core_witnesses(game: tuple[F, ...]):
    # A point game is not available here, so recover one feasible core point
    # by minimizing a zero quadratic from equal division.
    from scipy.optimize import linprog

    rows = []
    rhs = []
    for coalition in range(1, GRAND):
        rows.append(
            [
                -1.0 if coalition >> player & 1 else 0.0
                for player in range(N)
            ]
        )
        rhs.append(-float(game[coalition]))
    solved = linprog(
        np.zeros(N),
        A_ub=np.asarray(rows),
        b_ub=np.asarray(rhs),
        A_eq=np.ones((1, N)),
        b_eq=np.asarray([float(game[GRAND])]),
        bounds=[(None, None)] * N,
        method="highs",
    )
    if not solved.success:
        raise RuntimeError(solved.message)
    yield solved.x


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--points", type=int, default=8)
    parser.add_argument("--total", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    facets, _ = load_facets()
    tested = 0
    best = None
    cache: dict[tuple[F, ...], np.ndarray] = {}
    for trial in range(args.trials):
        points = {
            random_composition(N, args.total, rng) for _ in range(args.points)
        }
        if len(points) < 2:
            continue
        game = game_from_points(tuple(points))
        coalitions = list(range(1, GRAND))
        rng.shuffle(coalitions)
        for coalition in coalitions:
            cap = maximum_bump(game, coalition, facets)
            if cap <= 0:
                continue
            delta = cap * F(rng.randint(1, 9), 10)
            upper = list(game)
            upper[coalition] += delta
            upper_game = tuple(upper)
            if game not in cache:
                cache[game] = projected_shapley(game)
            if upper_game not in cache:
                cache[upper_game] = projected_shapley(upper_game)
            changes = cache[upper_game] - cache[game]
            minimum = min(
                float(changes[player])
                for player in range(N)
                if coalition >> player & 1
            )
            row = {
                "trial": trial,
                "coalition": coalition,
                "delta": str(delta),
                "minimum_protected_change": minimum,
                "changes": [float(value) for value in changes],
            }
            tested += 1
            if best is None or minimum < best[0]:
                best = (minimum, row, game, upper_game, tuple(points))
                print(json.dumps({"tested": tested, "best": row}), flush=True)
            if minimum < -1e-7:
                break
        if best is not None and best[0] < -1e-7:
            break
    payload = {
        "status": (
            "projected_shapley_counterexample_found"
            if best is not None and best[0] < -1e-7
            else "no_projected_shapley_counterexample_found"
        ),
        "tested_edges": tested,
        "configuration": vars(args) | {"output": str(args.output)},
        "best": None
        if best is None
        else best[1]
        | {
            "points": [[str(value) for value in point] for point in best[4]],
            "lower_game": [str(value) for value in best[2]],
            "upper_game": [str(value) for value in best[3]],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "tested_edges", "best")}, indent=2))
    return 1 if payload["status"] == "projected_shapley_counterexample_found" else 0


if __name__ == "__main__":
    raise SystemExit(main())
