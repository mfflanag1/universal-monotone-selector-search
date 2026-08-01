#!/usr/bin/env python3
"""Embed two differently oriented non-large exact games in a mixed four-cycle."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog

from n5_mixed_terminal_search import (
    GameCone,
    Topology,
    four_cycle_topology,
    lower_expectation_dual,
)


def nonlarge_game(permutation: tuple[int, ...]) -> np.ndarray:
    n = 5
    first = (3, 1, 1, 2, 2)
    second = (2, 3, 2, 1, 1)
    game = np.zeros(1 << n)
    for coalition in range(1 << n):
        game[coalition] = min(
            sum(
                first[permutation[player]]
                for player in range(n)
                if coalition >> player & 1
            ),
            sum(
                second[permutation[player]]
                for player in range(n)
                if coalition >> player & 1
            ),
        ) / 9.0
    return game


def shared_sink_topology() -> Topology:
    return Topology(
        5,
        ((0, 2, 1), (1, 2, 4), (0, 3, 2), (1, 4, 8)),
        (
            (2, 1, 0, 0, 0),
            (0, 0, 2, 1, 0),
            (-2, 0, -2, 0, 0),
            (0, -1, 0, 0, 0),
            (0, 0, 0, -1, 0),
        ),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", type=int, nargs=4, default=(2, 1, 2, 1))
    parser.add_argument(
        "--shape", choices=("four-cycle", "shared-sink"), default="four-cycle"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    topology = (
        shared_sink_topology()
        if args.shape == "shared-sink"
        else four_cycle_topology(tuple(args.weights), 5)
    )
    cone = GameCone(5, topology)
    first_game = nonlarge_game(tuple(range(5)))
    rows: list[dict[str, Any]] = []
    for permutation in itertools.permutations(range(5)):
        second_game = nonlarge_game(permutation)
        bounds: list[tuple[float | None, float | None]] = []
        for node in range(topology.node_count):
            for coalition in range(32):
                if node == 0:
                    value = float(first_game[coalition])
                    bounds.append((value, value))
                elif node == 1:
                    value = float(second_game[coalition])
                    bounds.append((value, value))
                else:
                    bounds.append((0.0, 1.0))
        result = linprog(
            np.zeros(cone.variable_count),
            A_ub=cone.a_ub,
            b_ub=cone.b_ub,
            A_eq=cone.a_eq,
            b_eq=cone.b_eq,
            bounds=bounds,
            method="highs",
        )
        if not result.success:
            continue
        games = [
            result.x[node * 32 : (node + 1) * 32]
            for node in range(topology.node_count)
        ]
        values = []
        for game, divergence in zip(
            games, topology.divergence, strict=True
        ):
            dual = lower_expectation_dual(game, divergence, 5)
            if dual is None:
                raise RuntimeError("lower-expectation dual failed")
            values.append(dual[0])
        rows.append(
            {
                "permutation": list(permutation),
                "objective": sum(values),
                "node_lower_expectations": values,
                "games": [[float(value) for value in game] for game in games],
            }
        )
    rows.sort(key=lambda row: float(row["objective"]), reverse=True)
    payload = {
        "status": "nonlarge_mixed_cycle_float",
        "weights": list(args.weights),
        "shape": args.shape,
        "feasible_orientations": len(rows),
        "positive_orientations": sum(
            float(row["objective"]) > 1e-9 for row in rows
        ),
        "best": rows[:20],
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
