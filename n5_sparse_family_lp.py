#!/usr/bin/env python3
"""Memory-safe sparse float LP for large exact-game families."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from family_search import Edge, Game, members


def sparse_family_result(
    games: Sequence[Game],
    edges: Sequence[Edge],
    n: int,
    method: str = "highs",
    dual_edge_weight_strategy: str | None = None,
    disp: bool = False,
    run_crossover: bool | None = None,
):
    allocation_variables = len(games) * n
    slack_index = allocation_variables
    variable_count = allocation_variables + 1
    grand = (1 << n) - 1

    row_indices: list[int] = []
    column_indices: list[int] = []
    coefficients: list[float] = []
    rhs: list[float] = []
    row = 0
    for game_index, game in enumerate(games):
        offset = game_index * n
        for coalition in range(1, grand):
            for player in members(coalition, n):
                row_indices.append(row)
                column_indices.append(offset + player)
                coefficients.append(-1.0)
            rhs.append(-float(game[coalition]))
            row += 1
    for lower, upper, coalition in edges:
        for player in members(coalition, n):
            row_indices.extend((row, row, row))
            column_indices.extend(
                (
                    lower * n + player,
                    upper * n + player,
                    slack_index,
                )
            )
            coefficients.extend((1.0, -1.0, 1.0))
            rhs.append(0.0)
            row += 1
    inequalities = coo_matrix(
        (
            coefficients,
            (row_indices, column_indices),
        ),
        shape=(row, variable_count),
    ).tocsr()

    equality_rows = np.repeat(np.arange(len(games)), n)
    equality_columns = np.arange(allocation_variables)
    equalities = coo_matrix(
        (
            np.ones(allocation_variables),
            (equality_rows, equality_columns),
        ),
        shape=(len(games), variable_count),
    ).tocsr()
    equality_rhs = np.array(
        [float(game[grand]) for game in games]
    )
    objective = np.zeros(variable_count)
    objective[slack_index] = -1.0
    options = {"disp": disp}
    if dual_edge_weight_strategy is not None:
        options["simplex_dual_edge_weight_strategy"] = dual_edge_weight_strategy
    if run_crossover is not None:
        options["run_crossover"] = "on" if run_crossover else "off"
    return linprog(
        objective,
        A_ub=inequalities,
        b_ub=np.array(rhs),
        A_eq=equalities,
        b_eq=equality_rhs,
        bounds=[(None, None)] * variable_count,
        method=method,
        options=options,
    )


def sparse_family_slack_float(
    games: Sequence[Game],
    edges: Sequence[Edge],
    n: int,
    method: str = "highs",
    dual_edge_weight_strategy: str | None = None,
    disp: bool = False,
    run_crossover: bool | None = None,
) -> float | None:
    result = sparse_family_result(
        games,
        edges,
        n,
        method,
        dual_edge_weight_strategy,
        disp,
        run_crossover,
    )
    if not result.success:
        return None
    return float(result.x[len(games) * n])
