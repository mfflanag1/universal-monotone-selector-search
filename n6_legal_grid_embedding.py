#!/usr/bin/env python3
"""Find a shared six-player exact extension of a four-player family."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from n6_mixed_terminal_search import load_facets


F = Fraction
N = 6
OLD_GRAND = 15
GRAND = 63
HELPERS = 48


def sparse(rows: list[dict[int, float]], columns: int):
    row_index = []
    column_index = []
    values = []
    for row, coefficients in enumerate(rows):
        for column, value in coefficients.items():
            if value:
                row_index.append(row)
                column_index.append(column)
                values.append(value)
    return coo_matrix(
        (values, (row_index, column_index)),
        shape=(len(rows), columns),
    ).tocsr()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.source.read_text())
    old_games = [
        [F(value) for value in game]
        for game in payload["family"]["games"]
    ]
    masks = [mask for mask in range(64) if mask & HELPERS]
    column = {mask: index for index, mask in enumerate(masks)}
    rows: list[dict[int, float]] = []
    rhs: list[float] = []
    facets = load_facets()
    for game in old_games:
        for facet in facets:
            row = {
                column[mask]: -float(facet[mask])
                for mask in masks
                if facet[mask]
            }
            old_part = sum(
                F(facet[mask]) * game[mask]
                for mask in range(16)
            )
            rows.append(row)
            rhs.append(float(old_part))
        for lower in range(64):
            for player in range(N):
                if lower >> player & 1:
                    continue
                upper = lower | (1 << player)
                lower_shared = bool(lower & HELPERS)
                upper_shared = bool(upper & HELPERS)
                if not lower_shared and not upper_shared:
                    if game[lower] > game[upper]:
                        raise RuntimeError("source game is nonmonotone")
                    continue
                row = {}
                bound = 0.0
                if lower_shared:
                    row[column[lower]] = 1.0
                else:
                    bound -= float(game[lower])
                if upper_shared:
                    row[column[upper]] = row.get(column[upper], 0.0) - 1.0
                else:
                    bound += float(game[upper])
                rows.append(row)
                rhs.append(bound)
    objective = np.zeros(len(masks))
    objective[column[GRAND]] = 1.0
    result = linprog(
        objective,
        A_ub=sparse(rows, len(masks)),
        b_ub=np.asarray(rhs),
        bounds=[(None, None)] * len(masks),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    shared = {
        mask: F(result.x[index]).limit_denominator(1_000_000)
        for mask, index in column.items()
    }
    games = []
    for old in old_games:
        game = [F(0)] * 64
        for mask in range(16):
            game[mask] = old[mask]
        for mask, value in shared.items():
            game[mask] = value
        games.append(game)
    edges = [
        [int(edge["lower"]), int(edge["upper"]), int(edge["coalition"])]
        for edge in payload["family"]["edges"]
    ]
    result_payload = {
        "status": "n6_shared_exact_extension_float",
        "source": str(args.source),
        "facet_count": len(facets),
        "game_count": len(games),
        "edge_count": len(edges),
        "grand_worth": str(shared[GRAND]),
        "shared_values": {str(mask): str(value) for mask, value in shared.items()},
        "games": [[str(value) for value in game] for game in games],
        "edges": edges,
    }
    args.output.write_text(json.dumps(result_payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result_payload["status"],
                "grand_worth": result_payload["grand_worth"],
                "constraints": len(rows),
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
