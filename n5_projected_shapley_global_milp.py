#!/usr/bin/env python3
"""Globally search exact n=5 edges for projected-Shapley monotonicity loss."""

from __future__ import annotations

import argparse
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
sys.path.insert(0, str(PROJECT_SRC))

from n5_facet_search import load_facets


N = 5
GRAND = (1 << N) - 1
COALITIONS = tuple(range(1, GRAND))
GAME_COUNT = GRAND + 1
DELTA = GAME_COUNT
LOWER_X = DELTA + 1
UPPER_X = LOWER_X + N
LOWER_LAMBDA = UPPER_X + N
UPPER_LAMBDA = LOWER_LAMBDA + len(COALITIONS)
LOWER_ACTIVE = UPPER_LAMBDA + len(COALITIONS)
UPPER_ACTIVE = LOWER_ACTIVE + len(COALITIONS)
VARIABLE_COUNT = UPPER_ACTIVE + len(COALITIONS)


def shapley_coefficients() -> np.ndarray:
    coefficients = np.zeros((N, GAME_COUNT))
    denominator = math.factorial(N)
    for player in range(N):
        for coalition in range(GAME_COUNT):
            size = coalition.bit_count()
            if coalition >> player & 1:
                coefficients[player, coalition] = (
                    math.factorial(size - 1)
                    * math.factorial(N - size)
                    / denominator
                )
            else:
                coefficients[player, coalition] = -(
                    math.factorial(size)
                    * math.factorial(N - size - 1)
                    / denominator
                )
    return coefficients


def centered_normal(coalition: int, player: int) -> float:
    return float(bool(coalition >> player & 1)) - coalition.bit_count() / N


def solve(
    coalition_size: int,
    minimum_bump: float,
    multiplier_bound: float,
    time_limit: float,
) -> dict[str, object]:
    changed = (1 << coalition_size) - 1
    protected_player = 0
    facets, _ = load_facets()
    shapley_matrix = shapley_coefficients()
    rows: list[dict[int, float]] = []
    lower: list[float] = []
    upper: list[float] = []

    def add(row, low=-np.inf, high=np.inf):
        rows.append({column: float(value) for column, value in row.items() if value})
        lower.append(float(low))
        upper.append(float(high))

    add({0: 1.0}, 0.0, 0.0)
    add({GRAND: 1.0}, 1.0, 1.0)
    for facet in facets:
        lower_row = {
            coalition: float(coefficient)
            for coalition, coefficient in enumerate(facet)
            if coefficient
        }
        add(lower_row, 0.0)
        upper_row = dict(lower_row)
        if facet[changed]:
            upper_row[DELTA] = float(facet[changed])
        add(upper_row, 0.0)
    for coalition in range(GAME_COUNT):
        for player in range(N):
            if coalition >> player & 1:
                continue
            successor = coalition | (1 << player)
            lower_row = {coalition: -1.0, successor: 1.0}
            add(lower_row, 0.0)
            upper_row = dict(lower_row)
            delta_coefficient = int(successor == changed) - int(coalition == changed)
            if delta_coefficient:
                upper_row[DELTA] = float(delta_coefficient)
            add(upper_row, 0.0)

    for side, x_start, lambda_start, active_start in (
        (0, LOWER_X, LOWER_LAMBDA, LOWER_ACTIVE),
        (1, UPPER_X, UPPER_LAMBDA, UPPER_ACTIVE),
    ):
        add({x_start + player: 1.0 for player in range(N)}, 1.0, 1.0)
        for player in range(N):
            stationarity = {x_start + player: 1.0}
            for worth, coefficient in enumerate(shapley_matrix[player]):
                if coefficient:
                    stationarity[worth] = stationarity.get(worth, 0.0) - coefficient
            if side:
                stationarity[DELTA] = -shapley_matrix[player, changed]
            for offset, coalition in enumerate(COALITIONS):
                normal = centered_normal(coalition, player)
                if normal:
                    stationarity[lambda_start + offset] = -normal
            add(stationarity, 0.0, 0.0)
        for offset, coalition in enumerate(COALITIONS):
            slack = {
                x_start + player: 1.0
                for player in range(N)
                if coalition >> player & 1
            }
            slack[coalition] = slack.get(coalition, 0.0) - 1.0
            if side and coalition == changed:
                slack[DELTA] = -1.0
            add(slack, 0.0)
            slack[active_start + offset] = float(N)
            add(slack, high=float(N))
            add(
                {
                    lambda_start + offset: 1.0,
                    active_start + offset: -multiplier_bound,
                },
                high=0.0,
            )
        # Conic Caratheodory in the four-dimensional efficiency hyperplane:
        # every normal-cone vector has a representation using at most four
        # linearly independent active core normals.
        add(
            {
                active_start + offset: 1.0
                for offset in range(len(COALITIONS))
            },
            high=float(N - 1),
        )

    rr: list[int] = []
    cc: list[int] = []
    vv: list[float] = []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            rr.append(row_index)
            cc.append(column)
            vv.append(value)
    matrix = coo_matrix(
        (vv, (rr, cc)), shape=(len(rows), VARIABLE_COUNT)
    ).tocsr()
    variable_lower = np.full(VARIABLE_COUNT, -np.inf)
    variable_upper = np.full(VARIABLE_COUNT, np.inf)
    variable_lower[:GAME_COUNT] = 0.0
    variable_upper[:GAME_COUNT] = 1.0
    variable_lower[DELTA] = minimum_bump
    variable_upper[DELTA] = 1.0
    variable_lower[LOWER_X : LOWER_X + 2 * N] = 0.0
    variable_upper[LOWER_X : LOWER_X + 2 * N] = 1.0
    variable_lower[LOWER_LAMBDA : LOWER_LAMBDA + 2 * len(COALITIONS)] = 0.0
    variable_upper[LOWER_LAMBDA : LOWER_LAMBDA + 2 * len(COALITIONS)] = multiplier_bound
    variable_lower[LOWER_ACTIVE:] = 0.0
    variable_upper[LOWER_ACTIVE:] = 1.0
    integrality = np.zeros(VARIABLE_COUNT)
    integrality[LOWER_ACTIVE:] = 1
    objective = np.zeros(VARIABLE_COUNT)
    objective[UPPER_X + protected_player] = 1.0
    objective[LOWER_X + protected_player] = -1.0
    result = milp(
        objective,
        integrality=integrality,
        bounds=Bounds(variable_lower, variable_upper),
        constraints=LinearConstraint(
            matrix, np.asarray(lower), np.asarray(upper)
        ),
        options={"time_limit": time_limit, "mip_rel_gap": 0.0},
    )
    payload: dict[str, object] = {
        "status": "optimal" if result.success else "incomplete",
        "message": result.message,
        "coalition_size": coalition_size,
        "changed_coalition": changed,
        "protected_player": protected_player,
        "minimum_bump": minimum_bump,
        "multiplier_bound": multiplier_bound,
        "objective": None if result.fun is None else float(result.fun),
        "mip_gap": getattr(result, "mip_gap", None),
        "node_count": getattr(result, "mip_node_count", None),
    }
    if result.x is not None:
        payload.update(
            {
                "delta": float(result.x[DELTA]),
                "lower_game": [float(value) for value in result.x[:GAME_COUNT]],
                "upper_game": [
                    float(result.x[coalition])
                    + (float(result.x[DELTA]) if coalition == changed else 0.0)
                    for coalition in range(GAME_COUNT)
                ],
                "lower_projection": [
                    float(value) for value in result.x[LOWER_X : LOWER_X + N]
                ],
                "upper_projection": [
                    float(value) for value in result.x[UPPER_X : UPPER_X + N]
                ],
                "maximum_multiplier": float(
                    np.max(
                        result.x[
                            LOWER_LAMBDA : LOWER_LAMBDA
                            + 2 * len(COALITIONS)
                        ]
                    )
                ),
                "active_lower": [
                    coalition
                    for offset, coalition in enumerate(COALITIONS)
                    if result.x[LOWER_ACTIVE + offset] > 0.5
                ],
                "active_upper": [
                    coalition
                    for offset, coalition in enumerate(COALITIONS)
                    if result.x[UPPER_ACTIVE + offset] > 0.5
                ],
            }
        )
    return payload


def rationalize(values: list[float], denominator: int) -> list[str]:
    return [str(Fraction(value).limit_denominator(denominator)) for value in values]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coalition-size", type=int, choices=range(1, N), required=True)
    parser.add_argument("--minimum-bump", type=float, default=0.001)
    parser.add_argument("--multiplier-bound", type=float, default=20.0)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--rational-denominator", type=int, default=100000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = solve(
        args.coalition_size,
        args.minimum_bump,
        args.multiplier_bound,
        args.time_limit,
    )
    if "lower_game" in payload:
        payload["rationalized"] = {
            key: rationalize(payload[key], args.rational_denominator)
            for key in (
                "lower_game",
                "upper_game",
                "lower_projection",
                "upper_projection",
            )
        }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    objective = payload["objective"]
    return 1 if objective is not None and objective < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
