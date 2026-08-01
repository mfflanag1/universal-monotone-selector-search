#!/usr/bin/env python3
"""Sample extreme legal exact edges and test projected Shapley values."""

from __future__ import annotations

import argparse
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
sys.path.insert(0, str(PROJECT_SRC))

from n5_facet_search import load_facets

from shapley_projection_edge_search import projected_shapley


N = 5
GRAND = (1 << N) - 1
GAME_COUNT = GRAND + 1
DELTA = GAME_COUNT


def legal_edge_polytope(changed: int, minimum_bump: float):
    facets, _ = load_facets()
    rows = []
    rhs = []
    for facet in facets:
        lower = np.zeros(GAME_COUNT + 1)
        lower[:GAME_COUNT] = -np.asarray(facet, dtype=float)
        rows.append(lower)
        rhs.append(0.0)
        upper = lower.copy()
        upper[DELTA] = -float(facet[changed])
        rows.append(upper)
        rhs.append(0.0)
    for coalition in range(GAME_COUNT):
        for player in range(N):
            if coalition >> player & 1:
                continue
            successor = coalition | (1 << player)
            lower = np.zeros(GAME_COUNT + 1)
            lower[coalition] = 1.0
            lower[successor] = -1.0
            rows.append(lower)
            rhs.append(0.0)
            upper = lower.copy()
            upper[DELTA] = float(coalition == changed) - float(
                successor == changed
            )
            rows.append(upper)
            rhs.append(0.0)
    equality = np.zeros((2, GAME_COUNT + 1))
    equality[0, 0] = 1.0
    equality[1, GRAND] = 1.0
    return (
        np.asarray(rows),
        np.asarray(rhs),
        equality,
        np.asarray([0.0, 1.0]),
        [(0.0, 1.0)] * GAME_COUNT + [(minimum_bump, 1.0)],
    )


def as_game(values: np.ndarray, denominator: int) -> tuple[Fraction, ...]:
    return tuple(
        Fraction(float(value)).limit_denominator(denominator)
        for value in values
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--minimum-bump", type=float, default=0.0001)
    parser.add_argument("--denominator", type=int, default=1000000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    best = None
    tested = 0
    for coalition_size in range(1, N):
        changed = (1 << coalition_size) - 1
        a_ub, b_ub, a_eq, b_eq, bounds = legal_edge_polytope(
            changed, args.minimum_bump
        )
        for sample in range(args.samples):
            objective = np.asarray(
                [rng.gauss(0.0, 1.0) for _ in range(GAME_COUNT + 1)]
            )
            result = linprog(
                objective,
                A_ub=a_ub,
                b_ub=b_ub,
                A_eq=a_eq,
                b_eq=b_eq,
                bounds=bounds,
                method="highs-ds",
            )
            if not result.success:
                raise RuntimeError(result.message)
            lower = as_game(result.x[:GAME_COUNT], args.denominator)
            upper = list(lower)
            delta = Fraction(float(result.x[DELTA])).limit_denominator(
                args.denominator
            )
            upper[changed] += delta
            upper_game = tuple(upper)
            lower_value = projected_shapley(lower)
            upper_value = projected_shapley(upper_game)
            changes = upper_value - lower_value
            minimum = min(changes[:coalition_size])
            tested += 1
            if best is None or minimum < best[0]:
                row = {
                    "coalition_size": coalition_size,
                    "changed_coalition": changed,
                    "sample": sample,
                    "delta": str(delta),
                    "minimum_protected_change": float(minimum),
                    "changes": [float(value) for value in changes],
                    "lower_game": [str(value) for value in lower],
                    "upper_game": [str(value) for value in upper_game],
                    "lower_projection": [float(value) for value in lower_value],
                    "upper_projection": [float(value) for value in upper_value],
                }
                best = (float(minimum), row)
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
        "best": None if best is None else best[1],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if payload["status"].endswith("_found") and not payload["status"].startswith("no_") else 0


if __name__ == "__main__":
    raise SystemExit(main())
