#!/usr/bin/env python3
"""Globally test protected-support upper preservation on the n=6 exact cone."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from n6_mixed_terminal_search import GRAND, N, load_facets


class Model:
    def __init__(self) -> None:
        self.game_offset = 0
        self.x_offset = GRAND + 1
        self.y_offset = self.x_offset + N
        self.lambda_offset = self.y_offset + N
        self.active_offset = self.lambda_offset + GRAND - 1
        self.beta = self.active_offset + GRAND - 1
        self.variable_count = self.beta + 1
        self.rows: list[dict[int, float]] = []
        self.lower: list[float] = []
        self.upper: list[float] = []

    def add(self, row: dict[int, float], lower=-np.inf, upper=np.inf) -> None:
        self.rows.append({column: value for column, value in row.items() if value})
        self.lower.append(lower)
        self.upper.append(upper)

    def matrix(self):
        rr: list[int] = []
        cc: list[int] = []
        vv: list[float] = []
        for row_index, row in enumerate(self.rows):
            for column, value in row.items():
                rr.append(row_index)
                cc.append(column)
                vv.append(value)
        return coo_matrix(
            (vv, (rr, cc)), shape=(len(self.rows), self.variable_count)
        ).tocsr()


def solve(support: int, weights: tuple[int, ...], time_limit: float) -> dict:
    model = Model()
    facets = load_facets()
    model.add({0: 1.0}, 0.0, 0.0)
    model.add({GRAND: 1.0}, 1.0, 1.0)
    for facet in facets:
        model.add(
            {
                coalition: float(coefficient)
                for coalition, coefficient in enumerate(facet)
                if coefficient
            },
            0.0,
        )
    for coalition in range(GRAND + 1):
        for player in range(N):
            if coalition >> player & 1:
                continue
            model.add(
                {
                    coalition | (1 << player): 1.0,
                    coalition: -1.0,
                },
                0.0,
            )

    model.add(
        {model.x_offset + player: 1.0 for player in range(N)}, 1.0, 1.0
    )
    model.add(
        {model.y_offset + player: 1.0 for player in range(N)}, 1.0, 1.0
    )
    for coalition in range(1, GRAND):
        core_row = {coalition: -1.0}
        for player in range(N):
            if coalition >> player & 1:
                core_row[model.x_offset + player] = 1.0
        model.add(core_row, 0.0)

        if coalition & support != support:
            relaxed_row = {coalition: -1.0}
            for player in range(N):
                if coalition >> player & 1:
                    relaxed_row[model.y_offset + player] = 1.0
            model.add(relaxed_row, 0.0)

        multiplier = model.lambda_offset + coalition - 1
        active = model.active_offset + coalition - 1
        model.add({multiplier: 1.0, active: -100.0}, upper=0.0)
        # x(T)-v(T) <= |T| (1-active).
        complementarity = {coalition: -1.0, active: float(coalition.bit_count())}
        for player in range(N):
            if coalition >> player & 1:
                complementarity[model.x_offset + player] = 1.0
        model.add(complementarity, upper=float(coalition.bit_count()))

    for player in range(N):
        stationarity = {model.beta: 1.0}
        for coalition in range(1, GRAND):
            if coalition >> player & 1:
                stationarity[model.lambda_offset + coalition - 1] = -1.0
        model.add(stationarity, float(weights[player]), float(weights[player]))

    objective = np.zeros(model.variable_count)
    for player, weight in enumerate(weights):
        objective[model.x_offset + player] = float(weight)
        objective[model.y_offset + player] = -float(weight)
    integrality = np.zeros(model.variable_count, dtype=np.uint8)
    integrality[model.active_offset : model.beta] = 1
    lower = np.zeros(model.variable_count)
    upper = np.ones(model.variable_count)
    upper[model.lambda_offset : model.active_offset] = 100.0
    lower[model.beta] = -100.0
    upper[model.beta] = 100.0
    result = milp(
        objective,
        integrality=integrality,
        bounds=Bounds(lower, upper),
        constraints=LinearConstraint(
            model.matrix(), np.asarray(model.lower), np.asarray(model.upper)
        ),
        options={"time_limit": time_limit, "mip_rel_gap": 0.0},
    )
    gap = None if result.fun is None else -float(result.fun)
    payload = {
        "status": result.message,
        "success": bool(result.success),
        "support": support,
        "weights": list(weights),
        "facet_count": len(facets),
        "gap": gap,
        "mip_gap": getattr(result, "mip_gap", None),
        "mip_node_count": getattr(result, "mip_node_count", None),
    }
    if result.x is not None:
        payload.update(
            {
                "game": [float(value) for value in result.x[: GRAND + 1]],
                "core_maximizer": [
                    float(value)
                    for value in result.x[model.x_offset : model.y_offset]
                ],
                "relaxed_maximizer": [
                    float(value)
                    for value in result.x[model.y_offset : model.lambda_offset]
                ],
            }
        )
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--support", type=int, required=True)
    parser.add_argument("--weights", type=int, nargs=N, required=True)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 0 < args.support < GRAND:
        raise ValueError("support must be nonempty and proper")
    if any(
        (args.weights[player] > 0) != bool(args.support >> player & 1)
        for player in range(N)
    ):
        raise ValueError("weights must be positive exactly on the support")
    payload = solve(args.support, tuple(args.weights), args.time_limit)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key != "game"}, indent=2))
    return 1 if payload.get("gap", 0.0) and payload["gap"] > 1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
