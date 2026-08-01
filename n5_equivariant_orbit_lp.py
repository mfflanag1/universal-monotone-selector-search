#!/usr/bin/env python3
"""Test exact-game families under the WLOG equivariance restriction."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_intrinsic_permutation_union import permute_game, permute_mask


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    n = int(payload["n"])
    grand = (1 << n) - 1
    if payload.get("family"):
        base_games = [
            tuple(F(value) for value in game)
            for game in payload["family"]["games"]
        ]
        base_edges = [
            (int(edge["lower"]), int(edge["upper"]), int(edge["coalition"]))
            for edge in payload["family"]["edges"]
        ]
    elif "base_v_by_mask" in payload:
        base_games = [
            tuple(F(value) for value in payload[key])
            for key in ("base_v_by_mask", "raised_v_by_mask")
        ]
        base_edges = [(0, 1, int(payload["bump"]["coalition_mask"]))]
    elif "base_worth" in payload:
        base_games = [
            tuple(F(str(value)) for value in payload[key])
            for key in ("base_worth", "raised_worth")
        ]
        base_edges = [(0, 1, int(payload["coalition"]))]
    else:
        raise ValueError("unrecognized family or two-game witness format")
    permutations = list(itertools.permutations(range(n)))
    games: list[tuple[Fraction, ...]] = []
    game_index: dict[tuple[Fraction, ...], int] = {}
    images: dict[tuple[int, tuple[int, ...]], int] = {}
    for base_index, game in enumerate(base_games):
        for permutation in permutations:
            image = permute_game(game, permutation)
            if image not in game_index:
                game_index[image] = len(games)
                games.append(image)
            images[base_index, permutation] = game_index[image]
    edges = set()
    for permutation in permutations:
        for lower, upper, coalition in base_edges:
            edges.add(
                (
                    images[lower, permutation],
                    images[upper, permutation],
                    permute_mask(coalition, permutation),
                )
            )
    allocation_count = len(games) * n
    slack_column = allocation_count
    variable_count = allocation_count + 1
    ub_rows: list[int] = []
    ub_columns: list[int] = []
    ub_values: list[float] = []
    ub_rhs: list[float] = []
    row = 0
    for index, game in enumerate(games):
        for coalition in range(1, grand):
            for player in range(n):
                if coalition >> player & 1:
                    ub_rows.append(row)
                    ub_columns.append(index * n + player)
                    ub_values.append(-1.0)
            ub_rhs.append(-float(game[coalition]))
            row += 1
    for lower, upper, coalition in sorted(edges):
        for player in range(n):
            if coalition >> player & 1:
                ub_rows.extend((row, row, row))
                ub_columns.extend(
                    (lower * n + player, upper * n + player, slack_column)
                )
                ub_values.extend((1.0, -1.0, 1.0))
                ub_rhs.append(0.0)
                row += 1
    a_ub = coo_matrix(
        (ub_values, (ub_rows, ub_columns)),
        shape=(row, variable_count),
    ).tocsr()
    eq_rows: list[int] = []
    eq_columns: list[int] = []
    eq_values: list[float] = []
    eq_rhs: list[float] = []
    row = 0
    for index, game in enumerate(games):
        for player in range(n):
            eq_rows.append(row)
            eq_columns.append(index * n + player)
            eq_values.append(1.0)
        eq_rhs.append(float(game[grand]))
        row += 1
    identity = tuple(range(n))
    for base_index in range(len(base_games)):
        base_image = images[base_index, identity]
        for permutation in permutations:
            image = images[base_index, permutation]
            for player in range(n):
                eq_rows.extend((row, row))
                eq_columns.extend(
                    (image * n + permutation[player], base_image * n + player)
                )
                eq_values.extend((1.0, -1.0))
                eq_rhs.append(0.0)
                row += 1
    a_eq = coo_matrix(
        (eq_values, (eq_rows, eq_columns)),
        shape=(row, variable_count),
    ).tocsr()
    objective = np.zeros(variable_count)
    objective[slack_column] = -1.0
    result = linprog(
        objective,
        A_ub=a_ub,
        b_ub=np.asarray(ub_rhs),
        A_eq=a_eq,
        b_eq=np.asarray(eq_rhs),
        bounds=[(None, None)] * variable_count,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    margin = float(result.x[slack_column])
    output = {
        "status": (
            "equivariant_obstruction_found"
            if margin < -1e-8
            else "equivariant_family_compatible"
        ),
        "source": str(args.input),
        "base_game_count": len(base_games),
        "base_edge_count": len(base_edges),
        "orbit_game_count": len(games),
        "orbit_edge_count": len(edges),
        "equivariance_equality_count": len(base_games) * len(permutations) * n,
        "max_min_margin_float": margin,
        "message": result.message,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
