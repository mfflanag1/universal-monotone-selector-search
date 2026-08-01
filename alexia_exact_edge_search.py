#!/usr/bin/env python3
"""Test the Alexia core selector on random legal exact-game bumps."""

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
from scipy.optimize import linprog

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_facet_search import load_facets

from family_search import game_from_points, random_composition


F = Fraction
N = 5
GRAND = (1 << N) - 1


def core_system(game: tuple[F, ...]):
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
    return np.asarray(rows), np.asarray(rhs)


def lexinal(game: tuple[F, ...], ordering: tuple[int, ...]) -> np.ndarray:
    a_ub, b_ub = core_system(game)
    equality_rows = [np.ones(N)]
    equality_rhs = [float(game[GRAND])]
    for player in ordering[:-1]:
        objective = np.zeros(N)
        objective[player] = -1.0
        result = linprog(
            objective,
            A_ub=a_ub,
            b_ub=b_ub,
            A_eq=np.asarray(equality_rows),
            b_eq=np.asarray(equality_rhs),
            bounds=[(None, None)] * N,
            method="highs",
        )
        if not result.success:
            raise RuntimeError(result.message)
        optimum = float(result.x[player])
        row = np.zeros(N)
        row[player] = 1.0
        equality_rows.append(row)
        equality_rhs.append(optimum)
    result = linprog(
        np.zeros(N),
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=np.asarray(equality_rows),
        b_eq=np.asarray(equality_rhs),
        bounds=[(None, None)] * N,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result.x


def alexia(game: tuple[F, ...]) -> np.ndarray:
    total = np.zeros(N)
    for ordering in itertools.permutations(range(N)):
        total += lexinal(game, ordering)
    return total / float(math.factorial(N))


def maximum_bump(
    game: tuple[F, ...], coalition: int, facets: list[tuple[int, ...]]
) -> F:
    bounds: list[F] = []
    for facet in facets:
        coefficient = facet[coalition]
        if coefficient < 0:
            slack = sum(
                F(facet[mask]) * game[mask]
                for mask in range(GRAND + 1)
            )
            bounds.append(slack / F(-coefficient))
    for player in range(N):
        if coalition >> player & 1:
            continue
        bounds.append(game[coalition | (1 << player)] - game[coalition])
    return min(bounds) if bounds else F(0)


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
    alexia_cache: dict[tuple[F, ...], np.ndarray] = {}
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
            if game not in alexia_cache:
                alexia_cache[game] = alexia(game)
            if upper_game not in alexia_cache:
                alexia_cache[upper_game] = alexia(upper_game)
            lower_value = alexia_cache[game]
            upper_value = alexia_cache[upper_game]
            changes = upper_value - lower_value
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
            if minimum < -1e-8:
                break
        if best is not None and best[0] < -1e-8:
            break
    payload = {
        "status": (
            "alexia_counterexample_found"
            if best is not None and best[0] < -1e-8
            else "no_alexia_counterexample_found"
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
    return 1 if payload["status"] == "alexia_counterexample_found" else 0


if __name__ == "__main__":
    raise SystemExit(main())
