#!/usr/bin/env python3
"""Solve a player-permutation-closed five-player family in symmetry quotient."""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from fractions import Fraction

import numpy as np
from ortools.sat.python import cp_model
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


N = 5
GRAND = (1 << N) - 1


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = bytearray(size)

    def find(self, value: int) -> int:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            parent = self.parent[value]
            self.parent[value] = root
            value = parent
        return root

    def union(self, left: int, right: int) -> None:
        left = self.find(left)
        right = self.find(right)
        if left == right:
            return
        if self.rank[left] < self.rank[right]:
            left, right = right, left
        self.parent[right] = left
        if self.rank[left] == self.rank[right]:
            self.rank[left] += 1


def adjacent_swap(player: int) -> tuple[list[int], list[int]]:
    player_map = list(range(N))
    player_map[player], player_map[player + 1] = player + 1, player
    coalition_map = []
    for coalition in range(GRAND + 1):
        image = 0
        for old_player in range(N):
            if coalition >> old_player & 1:
                image |= 1 << player_map[old_player]
        coalition_map.append(image)
    return player_map, coalition_map


def permute_game(game: bytes, coalition_map: list[int]) -> bytes:
    image = bytearray(GRAND + 1)
    for coalition, value in enumerate(game):
        image[coalition_map[coalition]] = value
    return bytes(image)


def quotient_variables(games: Sequence[bytes]) -> tuple[list[int], int]:
    index = {game: position for position, game in enumerate(games)}
    union_find = UnionFind(len(games) * N)
    for player in range(N - 1):
        player_map, coalition_map = adjacent_swap(player)
        for game_index, game in enumerate(games):
            image = index.get(permute_game(game, coalition_map))
            if image is None:
                raise ValueError("family is not closed under player permutations")
            for old_player in range(N):
                union_find.union(
                    game_index * N + old_player,
                    image * N + player_map[old_player],
                )
    roots: dict[int, int] = {}
    variables = []
    for raw in range(len(games) * N):
        root = union_find.find(raw)
        variables.append(roots.setdefault(root, len(roots)))
    return variables, len(roots)


def symmetric_family_result(
    games: Sequence[bytes],
    edges: Sequence[tuple[int, int, int]],
    method: str = "highs-ipm",
):
    variables, allocation_variable_count = quotient_variables(games)
    slack = allocation_variable_count

    core_rows: dict[tuple[tuple[int, int], ...], int] = {}
    for game_index, game in enumerate(games):
        offset = game_index * N
        for coalition in range(1, GRAND):
            counts = Counter(
                variables[offset + player]
                for player in range(N)
                if coalition >> player & 1
            )
            key = tuple(sorted(counts.items()))
            core_rows[key] = max(core_rows.get(key, game[coalition]), game[coalition])

    monotonicity_rows = set()
    for lower, upper, coalition in edges:
        for player in range(N):
            if coalition >> player & 1:
                monotonicity_rows.add(
                    (variables[lower * N + player], variables[upper * N + player])
                )

    row_indices: list[int] = []
    column_indices: list[int] = []
    coefficients: list[float] = []
    rhs: list[float] = []
    row = 0
    for key, worth in core_rows.items():
        for variable, count in key:
            row_indices.append(row)
            column_indices.append(variable)
            coefficients.append(-float(count))
        rhs.append(-float(worth))
        row += 1
    for lower, upper in monotonicity_rows:
        row_indices.extend((row, row, row))
        column_indices.extend((lower, upper, slack))
        coefficients.extend((1.0, -1.0, 1.0))
        rhs.append(0.0)
        row += 1
    inequalities = coo_matrix(
        (coefficients, (row_indices, column_indices)),
        shape=(row, allocation_variable_count + 1),
    ).tocsr()

    equality_rows: dict[tuple[tuple[int, int], ...], int] = {}
    for game_index, game in enumerate(games):
        counts = Counter(variables[game_index * N + player] for player in range(N))
        key = tuple(sorted(counts.items()))
        existing = equality_rows.get(key)
        if existing is not None and existing != game[GRAND]:
            raise ValueError("symmetry quotient merged inconsistent efficiency rows")
        equality_rows[key] = game[GRAND]
    eq_row_indices: list[int] = []
    eq_columns: list[int] = []
    eq_coefficients: list[float] = []
    equality_rhs = []
    for eq_row, (key, grand_worth) in enumerate(equality_rows.items()):
        for variable, count in key:
            eq_row_indices.append(eq_row)
            eq_columns.append(variable)
            eq_coefficients.append(float(count))
        equality_rhs.append(float(grand_worth))
    equalities = coo_matrix(
        (eq_coefficients, (eq_row_indices, eq_columns)),
        shape=(len(equality_rows), allocation_variable_count + 1),
    ).tocsr()

    objective = np.zeros(allocation_variable_count + 1)
    objective[slack] = -1.0
    result = linprog(
        objective,
        A_ub=inequalities,
        b_ub=np.asarray(rhs),
        A_eq=equalities,
        b_eq=np.asarray(equality_rhs),
        bounds=[(None, None)] * (allocation_variable_count + 1),
        method=method,
    )
    result.quotient_counts = {
        "allocation_variables": allocation_variable_count,
        "core_rows": len(core_rows),
        "monotonicity_rows": len(monotonicity_rows),
        "efficiency_rows": len(equality_rows),
    }
    result.quotient_variables = variables
    return result


def verify_rational_result(
    games: Sequence[bytes],
    edges: Sequence[tuple[int, int, int]],
    result,
    max_denominator: int = 1_000_000,
) -> dict[str, object]:
    values = [
        Fraction(float(value)).limit_denominator(max_denominator)
        for value in result.x
    ]
    allocation = values[:-1]
    margin = values[-1]
    variables = result.quotient_variables
    minimum_core_slack = None
    minimum_monotonicity_slack = None
    for game_index, game in enumerate(games):
        payoff = [allocation[variables[game_index * N + player]] for player in range(N)]
        if sum(payoff) != game[GRAND]:
            return {"verified": False, "reason": "efficiency", "game": game_index}
        for coalition in range(1, GRAND):
            slack = sum(
                payoff[player]
                for player in range(N)
                if coalition >> player & 1
            ) - game[coalition]
            minimum_core_slack = slack if minimum_core_slack is None else min(minimum_core_slack, slack)
            if slack < 0:
                return {
                    "verified": False,
                    "reason": "core",
                    "game": game_index,
                    "coalition": coalition,
                    "slack": str(slack),
                }
    for lower, upper, coalition in edges:
        for player in range(N):
            if not coalition >> player & 1:
                continue
            slack = (
                allocation[variables[upper * N + player]]
                - allocation[variables[lower * N + player]]
                - margin
            )
            minimum_monotonicity_slack = (
                slack
                if minimum_monotonicity_slack is None
                else min(minimum_monotonicity_slack, slack)
            )
            if slack < 0:
                return {
                    "verified": False,
                    "reason": "monotonicity",
                    "edge": [lower, upper, coalition],
                    "player": player,
                    "slack": str(slack),
                }
    return {
        "verified": True,
        "margin": str(margin),
        "minimum_core_slack": str(minimum_core_slack),
        "minimum_monotonicity_slack": str(minimum_monotonicity_slack),
        "nonzero_allocation_values": sum(value != 0 for value in allocation),
    }


def integer_grid_certificate(
    games: Sequence[bytes],
    edges: Sequence[tuple[int, int, int]],
    variables: Sequence[int],
    denominator: int = 2700,
    seconds: float = 120.0,
) -> dict[str, object]:
    variable_count = 1 + max(variables)
    maximum_grand = max(game[GRAND] for game in games)
    core_rows: dict[tuple[tuple[int, int], ...], int] = {}
    efficiency_rows: dict[tuple[tuple[int, int], ...], int] = {}
    for game_index, game in enumerate(games):
        game_variables = variables[game_index * N : (game_index + 1) * N]
        for coalition in range(1, GRAND):
            counts = Counter(
                game_variables[player]
                for player in range(N)
                if coalition >> player & 1
            )
            key = tuple(sorted(counts.items()))
            core_rows[key] = max(core_rows.get(key, game[coalition]), game[coalition])
        counts = Counter(game_variables)
        key = tuple(sorted(counts.items()))
        previous = efficiency_rows.get(key)
        if previous is not None and previous != game[GRAND]:
            raise ValueError("inconsistent efficiency row")
        efficiency_rows[key] = game[GRAND]
    monotonicity_rows = set()
    for lower, upper, coalition in edges:
        for player in range(N):
            if coalition >> player & 1:
                monotonicity_rows.add(
                    (variables[lower * N + player], variables[upper * N + player])
                )

    model = cp_model.CpModel()
    numerators = [
        model.new_int_var(0, denominator * maximum_grand, f"x_{index}")
        for index in range(variable_count)
    ]
    for key, worth in core_rows.items():
        model.add(sum(count * numerators[index] for index, count in key) >= denominator * worth)
    for key, grand_worth in efficiency_rows.items():
        model.add(sum(count * numerators[index] for index, count in key) == denominator * grand_worth)
    for lower, upper in monotonicity_rows:
        model.add(numerators[upper] - numerators[lower] >= 1)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 8
    status = solver.solve(model)
    payload: dict[str, object] = {
        "status": solver.status_name(status),
        "denominator": denominator,
        "margin": f"1/{denominator}",
        "variable_count": variable_count,
        "core_row_count": len(core_rows),
        "efficiency_row_count": len(efficiency_rows),
        "monotonicity_row_count": len(monotonicity_rows),
        "wall_time": solver.wall_time,
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        values = [solver.value(variable) for variable in numerators]
        payload["allocation_numerators"] = values
        payload["verified"] = all(
            sum(count * values[index] for index, count in key) >= denominator * worth
            for key, worth in core_rows.items()
        ) and all(
            sum(count * values[index] for index, count in key) == denominator * grand_worth
            for key, grand_worth in efficiency_rows.items()
        ) and all(values[upper] - values[lower] >= 1 for lower, upper in monotonicity_rows)
    else:
        payload["verified"] = False
    return payload
