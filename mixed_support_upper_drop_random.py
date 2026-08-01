#!/usr/bin/env python3
"""Randomly test protected-support upper preservation for exact games."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def optimum(
    game: np.ndarray,
    objective: np.ndarray,
    support: int | None,
) -> float:
    n = len(objective)
    grand = (1 << n) - 1
    rows = []
    rhs = []
    for coalition in range(1, grand):
        if support is not None and coalition & support == support:
            continue
        row = np.zeros(n)
        for player in range(n):
            if coalition >> player & 1:
                row[player] = -1.0
        rows.append(row)
        rhs.append(-game[coalition])
    result = linprog(
        -objective,
        A_ub=np.asarray(rows),
        b_ub=np.asarray(rhs),
        A_eq=np.ones((1, n)),
        b_eq=np.asarray([game[grand]]),
        bounds=[(None, None)] * n,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return -float(result.fun)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--cases", type=int, default=1000)
    parser.add_argument("--tests-per-case", type=int, default=100)
    parser.add_argument("--point-count", type=int, default=10)
    parser.add_argument("--total", type=int, default=100)
    parser.add_argument("--weight-bound", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    n = args.n
    grand = (1 << n) - 1
    rng = random.Random(args.seed)
    best: dict[str, object] | None = None
    tested = 0
    for case in range(args.cases):
        points = []
        for _ in range(args.point_count):
            cuts = sorted(rng.sample(range(1, args.total), n - 1))
            point = np.diff(np.asarray((0, *cuts, args.total), dtype=float))
            rng.shuffle(point)
            points.append(point)
        game = np.asarray(
            [
                min(
                    sum(
                        point[player]
                        for player in range(n)
                        if coalition >> player & 1
                    )
                    for point in points
                )
                for coalition in range(grand + 1)
            ]
        )
        for _ in range(args.tests_per_case):
            support_players = rng.sample(range(n), rng.randint(2, n - 1))
            support = sum(1 << player for player in support_players)
            objective = np.zeros(n)
            for player in support_players:
                objective[player] = rng.randint(1, args.weight_bound)
            if len(set(objective[player] for player in support_players)) == 1:
                continue
            original = optimum(game, objective, None)
            relaxed = optimum(game, objective, support)
            row = {
                "case": case,
                "support": support,
                "weights": [int(value) for value in objective],
                "original_upper": original,
                "relaxed_upper": relaxed,
                "gap": relaxed - original,
                "source_points": [
                    [int(value) for value in point] for point in points
                ],
            }
            tested += 1
            if best is None or float(row["gap"]) > float(best["gap"]):
                best = row
            if float(row["gap"]) > 1e-8:
                break
        if best is not None and float(best["gap"]) > 1e-8:
            break
        if (case + 1) % 100 == 0:
            print(
                json.dumps(
                    {
                        "cases": case + 1,
                        "tested": tested,
                        "best_gap": None if best is None else best["gap"],
                    }
                ),
                flush=True,
            )
    payload = {
        "status": (
            "mixed_upper_drop_found"
            if best is not None and float(best["gap"]) > 1e-8
            else "no_mixed_upper_drop_found"
        ),
        "n": n,
        "tested": tested,
        "best": best,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if payload["status"] == "mixed_upper_drop_found" else 0


if __name__ == "__main__":
    raise SystemExit(main())
