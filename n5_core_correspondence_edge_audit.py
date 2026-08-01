#!/usr/bin/env python3
"""Search exact bump pairs where a lower core point has no protected extension."""

from __future__ import annotations

import argparse
import json
import random
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


N = 5
GRAND = 31


def core_rows(game: list[float]):
    rows = []
    rhs = []
    for coalition in range(1, GRAND):
        rows.append(
            [-float(bool(coalition >> player & 1)) for player in range(N)]
        )
        rhs.append(-float(Fraction(game[coalition])))
    return rows, rhs


def optimize_core(game: list[float], objective: np.ndarray):
    rows, rhs = core_rows(game)
    return linprog(
        objective,
        A_ub=np.asarray(rows),
        b_ub=np.asarray(rhs),
        A_eq=np.ones((1, N)),
        b_eq=[float(Fraction(game[GRAND]))],
        bounds=[(None, None)] * N,
        method="highs",
    )


def extendable(game: list[float], point: np.ndarray, protected: int) -> bool:
    rows, rhs = core_rows(game)
    for player in range(N):
        if protected >> player & 1:
            row = [0.0] * N
            row[player] = -1.0
            rows.append(row)
            rhs.append(-float(point[player]))
    result = linprog(
        np.zeros(N),
        A_ub=np.asarray(rows),
        b_ub=np.asarray(rhs),
        A_eq=np.ones((1, N)),
        b_eq=[float(Fraction(game[GRAND]))],
        bounds=[(None, None)] * N,
        method="highs",
    )
    return result.success


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("family", type=Path)
    parser.add_argument("--edges", type=int, default=10000)
    parser.add_argument("--objectives", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = json.loads(args.family.read_text())
    raw = raw.get("family", raw)
    games = raw["games"]
    all_edges = [
        (
            int(edge["lower"]),
            int(edge["upper"]),
            int(edge["coalition"]),
        )
        if isinstance(edge, dict)
        else tuple(int(value) for value in edge)
        for edge in raw["edges"]
    ]
    rng = random.Random(args.seed)
    rng.shuffle(all_edges)
    tested_points = 0
    witness = None
    for edge_index, (lower, upper, coalition) in enumerate(
        all_edges[: args.edges]
    ):
        objectives = []
        for player in range(N):
            unit = np.zeros(N)
            unit[player] = -1.0
            objectives.append(unit)
            objectives.append(-unit)
        objectives.extend(
            np.asarray([rng.uniform(-1, 1) for _ in range(N)])
            for _ in range(args.objectives)
        )
        for objective in objectives:
            result = optimize_core(games[lower], objective)
            if not result.success:
                raise RuntimeError(result.message)
            tested_points += 1
            if not extendable(games[upper], result.x, coalition):
                witness = {
                    "edge": [lower, upper, coalition],
                    "lower_point": [float(value) for value in result.x],
                    "objective": [float(value) for value in objective],
                    "lower_game": games[lower],
                    "upper_game": games[upper],
                }
                break
        if witness is not None:
            break
        if edge_index % 1000 == 0:
            print(json.dumps({"edges": edge_index + 1, "points": tested_points}), flush=True)
    payload = {
        "status": "nonextendable_point_found" if witness else "no_nonextendable_point_found",
        "tested_edges": min(args.edges, len(all_edges)),
        "tested_points": tested_points,
        "witness": witness,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "tested_edges", "tested_points")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
