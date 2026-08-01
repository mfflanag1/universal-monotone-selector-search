#!/usr/bin/env python3
"""Report a floating common-core separation certificate for an archive."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from scipy.optimize import linprog


F = Fraction
N = 5
GRAND = (1 << N) - 1


def load_games(path: Path) -> list[list[float]]:
    payload = json.loads(path.read_text())
    if payload.get("family") is not None:
        return [
            [float(F(value)) for value in game]
            for game in payload["family"]["games"]
        ]
    record = payload.get("record", payload)
    return record["best_games_float"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", type=Path, nargs="+")
    args = parser.parse_args()

    for path in args.inputs:
        games = load_games(path)
        envelope = [
            max(game[coalition] for game in games)
            for coalition in range(1 << N)
        ]
        coalitions = list(range(1, GRAND))
        result = linprog(
            [1.0] * N,
            A_ub=[
                [
                    -1.0 if coalition >> player & 1 else 0.0
                    for player in range(N)
                ]
                for coalition in coalitions
            ],
            b_ub=[-envelope[coalition] for coalition in coalitions],
            bounds=[(None, None)] * N,
            method="highs",
        )
        if not result.success:
            raise RuntimeError(result.message)
        active = []
        for coalition, marginal in zip(
            coalitions, result.ineqlin.marginals, strict=True
        ):
            weight = -float(marginal)
            if weight <= 1e-9:
                continue
            maximum = envelope[coalition]
            nodes = [
                node
                for node, game in enumerate(games)
                if abs(game[coalition] - maximum) <= 1e-9
            ]
            active.append(
                {
                    "coalition": coalition,
                    "weight": weight,
                    "value": maximum,
                    "nodes": nodes,
                }
            )
        print(
            json.dumps(
                {
                    "file": str(path),
                    "gap": float(result.fun) - games[0][GRAND],
                    "allocation": list(result.x),
                    "active": active,
                }
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
