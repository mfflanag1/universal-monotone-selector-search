#!/usr/bin/env python3
"""Search for an equivariant-selector obstruction from a symmetric exact root."""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_facet_search import load_facets


N = 5
GRAND = (1 << N) - 1


def symmetric_exact_children(
    coalition: int, samples: int, rng: random.Random
) -> list[tuple[np.ndarray, float]]:
    facets, _ = load_facets()
    variable_count = GRAND + 2
    delta_column = GRAND + 1
    a_ub: list[np.ndarray] = []
    b_ub: list[float] = []
    a_eq: list[np.ndarray] = []
    b_eq: list[float] = []
    for facet in facets:
        root_row = np.zeros(variable_count)
        child_row = np.zeros(variable_count)
        for mask, coefficient in enumerate(facet):
            root_row[mask] = -float(coefficient)
            child_row[mask] = -float(coefficient)
        child_row[delta_column] = -float(facet[coalition])
        a_ub.extend((root_row, child_row))
        b_ub.extend((0.0, 0.0))
    for mask in range(GRAND + 1):
        for player in range(N):
            if mask >> player & 1:
                continue
            successor = mask | (1 << player)
            row = np.zeros(variable_count)
            row[mask] = 1.0
            row[successor] = -1.0
            a_ub.append(row)
            b_ub.append(0.0)
            if mask == coalition or successor == coalition:
                child_row = row.copy()
                child_row[delta_column] = (
                    1.0 if mask == coalition else -1.0
                )
                a_ub.append(child_row)
                b_ub.append(0.0)
    for mask in range(GRAND + 1):
        for other in range(mask + 1, GRAND + 1):
            if mask.bit_count() != other.bit_count():
                continue
            row = np.zeros(variable_count)
            row[mask] = 1.0
            row[other] = -1.0
            a_eq.append(row)
            b_eq.append(0.0)
    for mask, value in ((0, 0.0), (GRAND, 1.0)):
        row = np.zeros(variable_count)
        row[mask] = 1.0
        a_eq.append(row)
        b_eq.append(value)
    candidates = []
    objectives = []
    delta_objective = np.zeros(variable_count)
    delta_objective[delta_column] = -1.0
    objectives.append(delta_objective)
    for _ in range(samples):
        objective = np.zeros(variable_count)
        size_weights = [rng.uniform(-1.0, 1.0) for _ in range(N + 1)]
        for mask in range(GRAND + 1):
            objective[mask] = size_weights[mask.bit_count()]
        objective[delta_column] = rng.uniform(-4.0, -0.1)
        objectives.append(objective)
    for objective in objectives:
        result = linprog(
            objective,
            A_ub=np.asarray(a_ub),
            b_ub=np.asarray(b_ub),
            A_eq=np.asarray(a_eq),
            b_eq=np.asarray(b_eq),
            bounds=[(0.0, 1.0)] * (GRAND + 1) + [(1e-7, 1.0)],
            method="highs",
        )
        if result.success:
            candidates.append((result.x[: GRAND + 1], result.x[delta_column]))
    return candidates


def continuation_margin(root: np.ndarray, delta: float, coalition: int) -> float | None:
    # Maximize the common increment above the forced equal-root allocation for
    # members of the bumped coalition, over the child's core.
    t_column = N
    objective = np.zeros(N + 1)
    objective[t_column] = -1.0
    a_ub = []
    b_ub = []
    for mask in range(1, GRAND):
        row = np.zeros(N + 1)
        for player in range(N):
            if mask >> player & 1:
                row[player] = -1.0
        a_ub.append(row)
        b_ub.append(-float(root[mask] + (delta if mask == coalition else 0.0)))
    for player in range(N):
        if coalition >> player & 1:
            row = np.zeros(N + 1)
            row[player] = -1.0
            row[t_column] = 1.0
            a_ub.append(row)
            b_ub.append(-1.0 / N)
    equality = np.zeros((1, N + 1))
    equality[0, :N] = 1.0
    result = linprog(
        objective,
        A_ub=np.asarray(a_ub),
        b_ub=np.asarray(b_ub),
        A_eq=equality,
        b_eq=np.asarray([1.0]),
        bounds=[(None, None)] * (N + 1),
        method="highs",
    )
    if not result.success:
        return None
    return float(result.x[t_column])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    rows = []
    for size in range(1, N):
        coalition = (1 << size) - 1
        for root, delta in symmetric_exact_children(coalition, args.samples, rng):
            margin = continuation_margin(root, delta, coalition)
            if margin is None:
                continue
            rows.append(
                {
                    "coalition": coalition,
                    "coalition_size": size,
                    "delta": float(delta),
                    "continuation_margin": margin,
                    "root_by_size": [float(root[(1 << k) - 1]) for k in range(N + 1)],
                }
            )
    rows.sort(key=lambda row: float(row["continuation_margin"]))
    payload = {
        "status": (
            "equivariant_obstruction_found"
            if rows and float(rows[0]["continuation_margin"]) < -1e-8
            else "no_equivariant_obstruction_found"
        ),
        "tested": len(rows),
        "best": rows[:20],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
