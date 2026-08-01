#!/usr/bin/env python3
"""Maximize the mixed-divergence lower-expectation premium on exact games."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
sys.path.insert(0, str(PROJECT_SRC))

from n5_facet_search import load_facets
from n5_mixed_terminal_branch_search import extreme_representations


F = Fraction


def sparse_matrix(rows: list[dict[int, float]], columns: int):
    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[float] = []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                row_indices.append(row_index)
                column_indices.append(column)
                values.append(value)
    return coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(rows), columns),
    ).tocsr()


def choquet_representation(
    divergence: tuple[F, ...], n: int
) -> tuple[F, tuple[tuple[int, F], ...]]:
    levels = sorted(set(divergence))
    beta = levels[0]
    terms: list[tuple[int, F]] = []
    for lower, upper in zip(levels, levels[1:]):
        coalition = sum(
            1 << player
            for player, value in enumerate(divergence)
            if value >= upper
        )
        terms.append((coalition, upper - lower))
    return beta, tuple(terms)


def maximize_branch_premium(
    divergence: tuple[F, ...],
    representation: tuple[F, tuple[tuple[int, F], ...]],
    n: int,
    inequalities,
    rhs: np.ndarray,
    equalities: np.ndarray,
    equality_rhs: np.ndarray,
) -> tuple[float, np.ndarray] | None:
    grand = (1 << n) - 1
    beta, terms = representation
    choquet_beta, choquet_terms = choquet_representation(divergence, n)
    objective = np.zeros(grand + 1)
    objective[grand] += float(beta - choquet_beta)
    for coalition, weight in terms:
        objective[coalition] += float(weight)
    for coalition, weight in choquet_terms:
        objective[coalition] -= float(weight)
    result = linprog(
        -objective,
        A_ub=inequalities,
        b_ub=rhs,
        A_eq=equalities,
        b_eq=equality_rhs,
        bounds=[(None, None)] * (grand + 1),
        method="highs-ds",
        options={
            "dual_feasibility_tolerance": 1e-10,
            "primal_feasibility_tolerance": 1e-10,
        },
    )
    if not result.success:
        return None
    return -float(result.fun), result.x


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    n = 5
    grand = (1 << n) - 1
    facets, _ = load_facets()
    rows: list[dict[int, float]] = []
    rhs: list[float] = []
    for facet in facets:
        rows.append(
            {
                coalition: -float(coefficient)
                for coalition, coefficient in enumerate(facet)
                if coefficient
            }
        )
        rhs.append(0.0)
    for coalition in range(grand + 1):
        for player in range(n):
            if coalition >> player & 1:
                continue
            successor = coalition | (1 << player)
            rows.append({coalition: 1.0, successor: -1.0})
            rhs.append(0.0)
    equalities = np.zeros((2, grand + 1))
    equalities[0, 0] = 1.0
    equalities[1, grand] = 1.0
    equality_rhs = np.asarray([0.0, 1.0])
    inequality_matrix = sparse_matrix(rows, grand + 1)
    inequality_rhs = np.asarray(rhs)

    results: list[dict[str, Any]] = []
    for negative_count in range(1, n - 1):
        for zero_count in range(1, n - negative_count):
            positive_count = n - negative_count - zero_count
            divergence = (
                (F(-1),) * negative_count
                + (F(0),) * zero_count
                + (F(1),) * positive_count
            )
            representations = extreme_representations(list(divergence), n)
            best: dict[str, Any] | None = None
            for representation in representations:
                solved = maximize_branch_premium(
                    divergence,
                    representation,
                    n,
                    inequality_matrix,
                    inequality_rhs,
                    equalities,
                    equality_rhs,
                )
                if solved is None:
                    continue
                premium, game = solved
                if best is None or premium > float(best["premium"]):
                    beta, terms = representation
                    best = {
                        "premium": premium,
                        "beta": str(beta),
                        "terms": [
                            [coalition, str(weight)]
                            for coalition, weight in terms
                        ],
                        "game": [float(value) for value in game],
                    }
            if best is None:
                raise RuntimeError("no feasible lower-expectation branch")
            results.append(
                {
                    "level_multiplicities": [
                        negative_count,
                        zero_count,
                        positive_count,
                    ],
                    "divergence": [str(value) for value in divergence],
                    "representation_count": len(representations),
                    **best,
                }
            )
            print(
                (negative_count, zero_count, positive_count),
                len(representations),
                best["premium"],
            )
    payload = {
        "status": "five_player_three_level_premium_screen",
        "results": results,
    }
    if args.output:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
