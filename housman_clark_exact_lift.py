#!/usr/bin/env python3
"""Lift the Housman--Clark family to exact games using shared helper coalitions."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--grand", type=float, default=8.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.n < 5:
        raise ValueError("at least one helper player is required")

    raw = json.loads(args.source.read_text())
    old_games = [list(map(float, game)) for game in raw["family"]["games"]]
    node_count = len(old_games)
    n = args.n
    grand = (1 << n) - 1
    old_grand = 15
    helper_mask = grand ^ old_grand
    shared_masks = [mask for mask in range(1, grand + 1) if mask & helper_mask]
    shared_column = {mask: index for index, mask in enumerate(shared_masks)}
    witness_offset = len(shared_masks)
    proper = tuple(range(1, grand))

    def witness_column(node: int, tight: int, player: int) -> int:
        return witness_offset + ((node * len(proper) + tight - 1) * n + player)

    variable_count = witness_offset + node_count * len(proper) * n
    ub_rows: list[dict[int, float]] = []
    ub_rhs: list[float] = []
    eq_rows: list[dict[int, float]] = []
    eq_rhs: list[float] = []

    def game_term(node: int, mask: int) -> tuple[int | None, float]:
        if mask & helper_mask:
            return shared_column[mask], 0.0
        return None, old_games[node][mask]

    def add_game_term(row: dict[int, float], node: int, mask: int, coefficient: float) -> float:
        column, constant = game_term(node, mask)
        if column is None:
            return coefficient * constant
        row[column] = row.get(column, 0.0) + coefficient
        return 0.0

    # The helpers have zero singleton worth, the old four-player coalition and
    # the expanded grand coalition both have worth eight. Consequently every
    # core allocation gives every helper zero.
    for player in range(4, n):
        eq_rows.append({shared_column[1 << player]: 1.0})
        eq_rhs.append(0.0)
    eq_rows.append({shared_column[grand]: 1.0})
    eq_rhs.append(args.grand)

    for node in range(node_count):
        for lower in range(grand + 1):
            for player in range(n):
                if lower >> player & 1:
                    continue
                upper = lower | (1 << player)
                row: dict[int, float] = {}
                constant = add_game_term(row, node, lower, 1.0)
                constant += add_game_term(row, node, upper, -1.0)
                if row:
                    ub_rows.append(row)
                    ub_rhs.append(-constant)
                elif constant > 1e-9:
                    raise RuntimeError("the four-player source family is nonmonotone")

        for tight in proper:
            efficiency = {
                witness_column(node, tight, player): 1.0 for player in range(n)
            }
            eq_rows.append(efficiency)
            eq_rhs.append(args.grand)

            tightness: dict[int, float] = {}
            tight_constant = add_game_term(tightness, node, tight, 1.0)
            for player in range(n):
                if tight >> player & 1:
                    tightness[witness_column(node, tight, player)] = -1.0
            eq_rows.append(tightness)
            eq_rhs.append(-tight_constant)

            for coalition in proper:
                row: dict[int, float] = {}
                constant = add_game_term(row, node, coalition, 1.0)
                for player in range(n):
                    if coalition >> player & 1:
                        row[witness_column(node, tight, player)] = -1.0
                ub_rows.append(row)
                ub_rhs.append(-constant)

    def sparse(rows: list[dict[int, float]]):
        rr: list[int] = []
        cc: list[int] = []
        vv: list[float] = []
        for row_index, row in enumerate(rows):
            for column, value in row.items():
                if value:
                    rr.append(row_index)
                    cc.append(column)
                    vv.append(value)
        return coo_matrix((vv, (rr, cc)), shape=(len(rows), variable_count)).tocsr()

    result = linprog(
        np.zeros(variable_count),
        A_ub=sparse(ub_rows),
        b_ub=np.asarray(ub_rhs),
        A_eq=sparse(eq_rows),
        b_eq=np.asarray(eq_rhs),
        bounds=[(0.0, args.grand)] * variable_count,
        method="highs",
    )
    payload: dict[str, object] = {
        "status": "feasible" if result.success else "infeasible",
        "message": result.message,
        "n": n,
        "node_count": node_count,
        "shared_game_variables": len(shared_masks),
        "witness_variables": variable_count - witness_offset,
        "inequalities": len(ub_rows),
        "equalities": len(eq_rows),
    }
    if result.success:
        shared = {mask: result.x[column] for mask, column in shared_column.items()}
        games = []
        for node in range(node_count):
            game = []
            for mask in range(grand + 1):
                game.append(shared[mask] if mask & helper_mask else old_games[node][mask])
            games.append(game)
        payload.update(
            {
                "games": games,
                "edges": [
                    [int(edge["lower"]), int(edge["upper"]), int(edge["coalition"])]
                    for edge in raw["family"]["edges"]
                ],
                "shared_values": {str(mask): value for mask, value in shared.items()},
            }
        )
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key not in {"games", "shared_values"}}, indent=2))
    return 0 if result.success else 2


if __name__ == "__main__":
    raise SystemExit(main())
