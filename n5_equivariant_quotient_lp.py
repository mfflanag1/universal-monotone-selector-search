#!/usr/bin/env python3
"""Permutation-quotient LP for equivariant coalitionally monotone core selection."""

from __future__ import annotations

import itertools

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from n5_intrinsic_permutation_union import permute_game, permute_mask


def canonical_game(game, permutations):
    images = [(permute_game(game, permutation), permutation) for permutation in permutations]
    return min(images, key=lambda item: item[0])


def permutation_indices(n, permutations):
    size = 1 << n
    rows = []
    for permutation in permutations:
        inverse = [0] * size
        for coalition in range(size):
            inverse[permute_mask(coalition, permutation)] = coalition
        rows.append(tuple(inverse))
    return rows


def bulk_canonical_games(games, n, permutations):
    indices = permutation_indices(n, permutations)
    values = np.asarray([[float(value) for value in game] for game in games])
    best_values = values[:, indices[0]].copy()
    best_permutations = np.zeros(len(games), dtype=np.int16)
    for permutation_index, image_indices in enumerate(indices[1:], start=1):
        candidate = values[:, image_indices]
        undecided = np.ones(len(games), dtype=bool)
        less = np.zeros(len(games), dtype=bool)
        for column in range(values.shape[1]):
            less |= undecided & (candidate[:, column] < best_values[:, column])
            undecided &= candidate[:, column] == best_values[:, column]
            if not np.any(undecided):
                break
        if np.any(less):
            best_values[less] = candidate[less]
            best_permutations[less] = permutation_index
    canonical_games = []
    canonical_index = {}
    maps = []
    for game, permutation_index in zip(games, best_permutations, strict=True):
        canonical = tuple(game[index] for index in indices[int(permutation_index)])
        if canonical not in canonical_index:
            canonical_index[canonical] = len(canonical_games)
            canonical_games.append(canonical)
        maps.append(
            (canonical_index[canonical], permutations[int(permutation_index)])
        )
    canonical_values = np.asarray(
        [[float(value) for value in game] for game in canonical_games]
    )
    stabilizers = [[] for _ in canonical_games]
    for permutation, image_indices in zip(permutations, indices, strict=True):
        fixed = np.all(
            canonical_values[:, image_indices] == canonical_values,
            axis=1,
        )
        for game_index in np.flatnonzero(fixed):
            stabilizers[int(game_index)].append(permutation)
    return canonical_games, maps, stabilizers


def quotient_family(games, edges, n=5):
    permutations = list(itertools.permutations(range(n)))
    canonical_games, maps, stabilizers = bulk_canonical_games(
        games, n, permutations
    )
    quotient_edges = set()
    for lower, upper, coalition in edges:
        lower_orbit, lower_map = maps[lower]
        upper_orbit, upper_map = maps[upper]
        for player in range(n):
            if coalition >> player & 1:
                quotient_edges.add(
                    (
                        lower_orbit,
                        lower_map[player],
                        upper_orbit,
                        upper_map[player],
                    )
                )
    return canonical_games, quotient_edges, stabilizers


def solve_quotient(
    canonical_games,
    quotient_edges,
    stabilizers,
    n=5,
    method="highs-ipm",
    return_result=False,
):
    variable_count = len(canonical_games) * n
    slack = variable_count
    ub_rows = []
    ub_columns = []
    ub_values = []
    ub_rhs = []
    row = 0
    grand = (1 << n) - 1
    for game_index, game in enumerate(canonical_games):
        for coalition in range(1, grand):
            for player in range(n):
                if coalition >> player & 1:
                    ub_rows.append(row)
                    ub_columns.append(game_index * n + player)
                    ub_values.append(-1.0)
            ub_rhs.append(-float(game[coalition]))
            row += 1
    for lower_orbit, lower_player, upper_orbit, upper_player in sorted(quotient_edges):
        ub_rows.extend((row, row, row))
        ub_columns.extend(
            (
                lower_orbit * n + lower_player,
                upper_orbit * n + upper_player,
                slack,
            )
        )
        ub_values.extend((1.0, -1.0, 1.0))
        ub_rhs.append(0.0)
        row += 1
    inequalities = coo_matrix(
        (ub_values, (ub_rows, ub_columns)),
        shape=(row, variable_count + 1),
    ).tocsr()

    eq_rows = []
    eq_columns = []
    eq_values = []
    eq_rhs = []
    row = 0
    for game_index, game in enumerate(canonical_games):
        for player in range(n):
            eq_rows.append(row)
            eq_columns.append(game_index * n + player)
            eq_values.append(1.0)
        eq_rhs.append(float(game[grand]))
        row += 1
        for permutation in stabilizers[game_index]:
            for player in range(n):
                if permutation[player] == player:
                    continue
                eq_rows.extend((row, row))
                eq_columns.extend(
                    (
                        game_index * n + permutation[player],
                        game_index * n + player,
                    )
                )
                eq_values.extend((1.0, -1.0))
                eq_rhs.append(0.0)
                row += 1
    equalities = coo_matrix(
        (eq_values, (eq_rows, eq_columns)),
        shape=(row, variable_count + 1),
    ).tocsr()
    objective = np.zeros(variable_count + 1)
    objective[slack] = -1.0
    result = linprog(
        objective,
        A_ub=inequalities,
        b_ub=np.asarray(ub_rhs),
        A_eq=equalities,
        b_eq=np.asarray(eq_rhs),
        bounds=[(None, None)] * (variable_count + 1),
        method=method,
    )
    if not result.success:
        raise RuntimeError(result.message)
    margin = float(result.x[slack])
    return (margin, result) if return_result else margin


def equivariant_quotient_margin(games, edges, n=5, method="highs-ipm"):
    """Return the maximum common CM margin on the permutation closure of a family."""
    canonical_games, quotient_edges, stabilizers = quotient_family(
        games, edges, n
    )
    margin = solve_quotient(
        canonical_games, quotient_edges, stabilizers, n, method
    )
    return margin, len(canonical_games), len(quotient_edges)
