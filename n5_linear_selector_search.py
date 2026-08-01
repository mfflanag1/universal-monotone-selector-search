#!/usr/bin/env python3
"""Find a coalitionally monotone linear core selector on the n=5 exact cone."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from n5_facet_search import load_facets


F = Fraction
N = 5
GRAND = (1 << N) - 1
COALITIONS = tuple(range(1, GRAND + 1))
PROPER = tuple(range(1, GRAND))


def a_index(player: int, coalition: int) -> int:
    return player * GRAND + coalition - 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-denominator", type=int, default=10**7)
    parser.add_argument("--method", default="highs")
    parser.add_argument("--no-presolve", action="store_true")
    parser.add_argument("--include-monotone-game", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    facets, _ = load_facets()
    if args.include_monotone_game:
        for coalition in range(GRAND + 1):
            for player in range(N):
                if coalition >> player & 1:
                    continue
                successor = coalition | (1 << player)
                row = [0] * (GRAND + 1)
                if coalition:
                    row[coalition] = -1
                row[successor] = 1
                facets.append(tuple(row))
    a_count = N * GRAND
    lambda_offset = a_count

    def lambda_index(core_coalition: int, facet: int) -> int:
        return lambda_offset + (core_coalition - 1) * len(facets) + facet

    variable_count = a_count + len(PROPER) * len(facets)
    rows: list[dict[int, float]] = []
    rhs: list[float] = []

    # Efficiency as an identity in the game coordinates.
    for coalition in COALITIONS:
        rows.append(
            {a_index(player, coalition): 1.0 for player in range(N)}
        )
        rhs.append(1.0 if coalition == GRAND else 0.0)

    # Each core slack is a nonnegative combination of exact-cone facets.
    for core_coalition in PROPER:
        members = [
            player
            for player in range(N)
            if core_coalition >> player & 1
        ]
        for game_coalition in COALITIONS:
            row = {
                a_index(player, game_coalition): 1.0
                for player in members
            }
            for facet_index, facet in enumerate(facets):
                coefficient = facet[game_coalition]
                if coefficient:
                    row[lambda_index(core_coalition, facet_index)] = -float(
                        coefficient
                    )
            rows.append(row)
            rhs.append(1.0 if game_coalition == core_coalition else 0.0)

    row_index: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                row_index.append(index)
                columns.append(column)
                values.append(value)
    matrix = coo_matrix(
        (values, (row_index, columns)),
        shape=(len(rows), variable_count),
    ).tocsr()
    bounds = [(None, None)] * a_count + [(0.0, None)] * (
        variable_count - a_count
    )
    # Coalitional monotonicity fixes the signs of the relevant A entries.
    for coalition in PROPER:
        for player in range(N):
            if coalition >> player & 1:
                bounds[a_index(player, coalition)] = (0.0, None)
    result = linprog(
        np.zeros(variable_count),
        A_eq=matrix,
        b_eq=np.asarray(rhs),
        bounds=bounds,
        method=args.method,
        options={"presolve": not args.no_presolve},
    )
    payload: dict[str, object] = {
        "status": "feasible" if result.success else "infeasible",
        "message": result.message,
        "facet_count": len(facets),
        "variable_count": variable_count,
        "equality_count": len(rows),
    }
    if result.success:
        rational = [
            F(float(value)).limit_denominator(args.max_denominator)
            for value in result.x
        ]
        exact_equalities = True
        for row, expected in zip(rows, rhs, strict=True):
            actual = sum(
                F(str(value)) * rational[column]
                for column, value in row.items()
            )
            if actual != F(str(expected)):
                exact_equalities = False
                break
        exact_bounds = all(
            (lower is None or rational[index] >= F(str(lower)))
            and (upper is None or rational[index] <= F(str(upper)))
            for index, (lower, upper) in enumerate(bounds)
        )
        selector = [
            [str(rational[a_index(player, coalition)]) for coalition in COALITIONS]
            for player in range(N)
        ]
        payload.update(
            {
                "exact_reconstruction": exact_equalities and exact_bounds,
                "selector_coefficients": selector,
                "nonzero_selector_coefficients": sum(
                    value != 0 for value in rational[:a_count]
                ),
                "minimum_protected_coefficient": str(
                    min(
                        rational[a_index(player, coalition)]
                        for coalition in PROPER
                        for player in range(N)
                        if coalition >> player & 1
                    )
                ),
            }
        )
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0 if result.success else 1


if __name__ == "__main__":
    raise SystemExit(main())
