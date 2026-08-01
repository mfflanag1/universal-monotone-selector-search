#!/usr/bin/env python3
"""Audit a fixed-priority lexicographic core selector on an exact family."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def lexicographic_core(game, order: tuple[int, ...]) -> list[float]:
    n = len(order)
    grand = (1 << n) - 1
    rows = []
    rhs = []
    for coalition in range(1, grand):
        rows.append(
            [-float(bool(coalition >> player & 1)) for player in range(n)]
        )
        rhs.append(-float(Fraction(game[coalition])))
    equality_rows = [np.ones(n)]
    equality_rhs = [float(Fraction(game[grand]))]
    point = None
    for player in order:
        objective = np.zeros(n)
        objective[player] = -1.0
        result = linprog(
            objective,
            A_ub=np.asarray(rows),
            b_ub=np.asarray(rhs),
            A_eq=np.asarray(equality_rows),
            b_eq=np.asarray(equality_rhs),
            bounds=[(None, None)] * n,
            method="highs",
        )
        if not result.success:
            raise RuntimeError(result.message)
        point = result.x
        equality = np.zeros(n)
        equality[player] = 1.0
        equality_rows.append(equality)
        equality_rhs.append(float(point[player]))
    if point is None:
        raise RuntimeError("empty player order")
    return [float(value) for value in point]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("family", type=Path)
    parser.add_argument("--order", default="0,1,2,3,4")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = json.loads(args.family.read_text())
    raw = raw.get("family", raw)
    games = raw["games"]
    order = tuple(int(value) for value in args.order.split(","))
    allocations = [lexicographic_core(game, order) for game in games]
    violations = []
    for edge in raw["edges"]:
        if isinstance(edge, dict):
            lower, upper, coalition = (
                int(edge["lower"]),
                int(edge["upper"]),
                int(edge["coalition"]),
            )
        else:
            lower, upper, coalition = (int(value) for value in edge)
        for player in range(len(order)):
            difference = allocations[upper][player] - allocations[lower][player]
            if coalition >> player & 1 and difference < -1e-8:
                violations.append(
                    {
                        "edge": [lower, upper, coalition],
                        "player": player,
                        "difference": difference,
                        "lower_allocation": allocations[lower],
                        "upper_allocation": allocations[upper],
                    }
                )
                break
    payload = {
        "status": "violation_found" if violations else "clean",
        "order": list(order),
        "games": len(games),
        "edges": len(raw["edges"]),
        "violations": violations[:20],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
