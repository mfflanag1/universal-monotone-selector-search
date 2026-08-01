#!/usr/bin/env python3
"""Search finite exact-game families for monotone-selection obstructions.

For games ``v^0, ..., v^(m-1)`` and directed bump edges ``(k, l, S)``,
where ``v^l = v^k + d 1_S``, a coalitionally monotonic core selector would
have to choose allocations ``x^k`` satisfying

    x^k in core(v^k), and x^k_i <= x^l_i for every i in S.

These conditions form one linear feasibility problem.  An infeasible family
would disprove the existence of a coalitionally monotonic core selection on
any domain containing it.

The search starts from exact games generated as lower envelopes of common-sum
integer points.  It then grows a directed acyclic graph using single-coalition
bumps that preserve exactness.  This is deliberately stronger than adding
independent private points to one common upper game: that simpler construction
has a common core point and can never be an obstruction.

SciPy/HiGHS is used only for screening.  Any apparent infeasible family is
rechecked with the exact simplex solver and an exact Farkas certificate.
"""

from __future__ import annotations

import argparse
import json
import random
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence

from exact_lp import exact_linprog


F = Fraction
Game = tuple[F, ...]
Point = tuple[F, ...]
Edge = tuple[int, int, int]


def exact_dot(left: Sequence[F], right: Sequence[F]) -> F:
    """Return an exact dot product without multiplying structural zeros."""

    return sum(
        (
            coefficient * value
            for coefficient, value in zip(left, right, strict=True)
            if coefficient
        ),
        start=F(0),
    )


def scipy_linprog(*args: Any, **kwargs: Any) -> Any:
    """Import SciPy only on float-screening paths.

    Exact verification remains standalone and needs only the standard library.
    """

    from scipy.optimize import linprog

    return linprog(*args, **kwargs)


def members(mask: int, n: int) -> list[int]:
    return [player for player in range(n) if mask >> player & 1]


def coalition_sum(point: Sequence[F], mask: int) -> F:
    return sum(
        (point[player] for player in members(mask, len(point))),
        start=F(0),
    )


def random_composition(
    n: int, total: int, rng: random.Random
) -> Point:
    cuts = sorted(rng.randint(0, total) for _ in range(n - 1))
    values: list[int] = []
    previous = 0
    for cut in cuts:
        values.append(cut - previous)
        previous = cut
    values.append(total - previous)
    rng.shuffle(values)
    return tuple(F(value) for value in values)


def simplex_points(n: int, total: int) -> list[Point]:
    return [
        tuple(
            F(total) if player == support else F(0)
            for player in range(n)
        )
        for support in range(n)
    ]


def near_simplex_points(
    n: int, total: int, rng: random.Random
) -> list[Point]:
    noise_max = max(1, total // (2 * n))
    points: list[Point] = []
    for support in range(n):
        values = [
            rng.randint(0, noise_max) if player != support else 0
            for player in range(n)
        ]
        noise = sum(values)
        if noise > total:
            values = [0] * n
            noise = 0
        values[support] = total - noise
        points.append(tuple(F(value) for value in values))
    return points


def game_from_points(points: Sequence[Point]) -> Game:
    n = len(points[0])
    grand = (1 << n) - 1
    total = sum(points[0])
    assert all(sum(point) == total for point in points)
    worth = [F(0)] * (1 << n)
    for coalition in range(1, grand):
        worth[coalition] = min(
            coalition_sum(point, coalition) for point in points
        )
    worth[grand] = total
    return tuple(worth)


def bump_game(game: Game, coalition: int, delta: F) -> Game:
    bumped = list(game)
    bumped[coalition] += delta
    return tuple(bumped)


def is_monotone(game: Game, n: int) -> bool:
    grand = (1 << n) - 1
    return all(
        game[coalition] <= game[coalition | (1 << player)]
        for coalition in range(grand + 1)
        for player in range(n)
        if not coalition >> player & 1
    )


def is_superadditive(game: Game, n: int) -> bool:
    grand = (1 << n) - 1
    for left in range(grand + 1):
        remaining = grand ^ left
        right = remaining
        while True:
            if game[left | right] < game[left] + game[right]:
                return False
            if right == 0:
                break
            right = (right - 1) & remaining
    return True


def core_rows(
    game: Game, n: int
) -> tuple[list[list[float]], list[float]]:
    grand = (1 << n) - 1
    rows: list[list[float]] = []
    rhs: list[float] = []
    for coalition in range(1, grand):
        rows.append(
            [
                -1.0 if coalition >> player & 1 else 0.0
                for player in range(n)
            ]
        )
        rhs.append(-float(game[coalition]))
    return rows, rhs


def tight_core_point_float(
    game: Game, n: int, tight_coalition: int
) -> list[float] | None:
    rows, rhs = core_rows(game, n)
    equality_rows = [[1.0] * n]
    equality_rhs = [float(game[(1 << n) - 1])]
    if tight_coalition != (1 << n) - 1:
        equality_rows.append(
            [
                1.0 if tight_coalition >> player & 1 else 0.0
                for player in range(n)
            ]
        )
        equality_rhs.append(float(game[tight_coalition]))
    result = scipy_linprog(
        [0.0] * n,
        A_ub=rows,
        b_ub=rhs,
        A_eq=equality_rows,
        b_eq=equality_rhs,
        bounds=[(None, None)] * n,
        method="highs",
    )
    if not result.success:
        return None
    return [float(value) for value in result.x]


def is_exact_float(game: Game, n: int) -> bool:
    grand = (1 << n) - 1
    return all(
        tight_core_point_float(game, n, coalition) is not None
        for coalition in range(1, grand)
    )


def tight_core_point_exact(
    game: Game, n: int, tight_coalition: int
) -> Point | None:
    grand = (1 << n) - 1
    inequalities: list[list[F]] = []
    inequality_rhs: list[F] = []
    for coalition in range(1, grand):
        inequalities.append(
            [
                F(-1) if coalition >> player & 1 else F(0)
                for player in range(n)
            ]
        )
        inequality_rhs.append(-game[coalition])
    equalities = [[F(1)] * n]
    equality_rhs = [game[grand]]
    if tight_coalition != grand:
        equalities.append(
            [
                F(1) if tight_coalition >> player & 1 else F(0)
                for player in range(n)
            ]
        )
        equality_rhs.append(game[tight_coalition])
    status, _, point = exact_linprog(
        [F(0)] * n,
        inequalities,
        inequality_rhs,
        equalities,
        equality_rhs,
    )
    if status != "optimal" or point is None:
        return None
    return tuple(point)


def exactness_certificate(game: Game, n: int) -> list[Point] | None:
    grand = (1 << n) - 1
    witnesses: list[Point] = []
    for coalition in range(1, grand):
        point = tight_core_point_exact(game, n, coalition)
        if point is None:
            return None
        witnesses.append(point)
    return witnesses


def check_tight_core_point(
    game: Game, point: Point, n: int, tight_coalition: int
) -> bool:
    grand = (1 << n) - 1
    return (
        sum(point) == game[grand]
        and coalition_sum(point, tight_coalition)
        == game[tight_coalition]
        and all(
            coalition_sum(point, coalition) >= game[coalition]
            for coalition in range(1, grand)
        )
    )


def reconstructed_exactness_certificate(
    game: Game, n: int, max_denominator: int = 1_000_000
) -> list[Point] | None:
    """Reconstruct float LP vertices, then verify them exactly.

    This keeps bulk archival verification fast without trusting floating-point
    feasibility.  Failure to reconstruct is inconclusive and callers may fall
    back to the exact simplex.
    """

    grand = (1 << n) - 1
    witnesses: list[Point] = []
    for coalition in range(1, grand):
        float_point = tight_core_point_float(game, n, coalition)
        if float_point is None:
            return None
        point = tuple(
            F(value).limit_denominator(max_denominator)
            for value in float_point
        )
        if not check_tight_core_point(game, point, n, coalition):
            return None
        witnesses.append(point)
    return witnesses


def arm_pool(
    game: Game,
    n: int,
    bump_unit: F,
    bump_steps: int,
    coalitions_per_node: int,
    require_monotone: bool,
    require_superadditive: bool,
    rng: random.Random,
) -> list[tuple[int, F, Game]]:
    grand = (1 << n) - 1
    coalitions = list(range(1, grand))
    rng.shuffle(coalitions)
    if coalitions_per_node > 0:
        coalitions = coalitions[:coalitions_per_node]

    arms: list[tuple[int, F, Game]] = []
    for coalition in coalitions:
        best: tuple[int, F, Game] | None = None
        for step in range(1, bump_steps + 1):
            delta = bump_unit * step
            candidate = bump_game(game, coalition, delta)
            if require_monotone and not is_monotone(candidate, n):
                break
            if require_superadditive and not is_superadditive(candidate, n):
                break
            if not is_exact_float(candidate, n):
                break
            best = coalition, delta, candidate
        if best is not None:
            arms.append(best)
    rng.shuffle(arms)
    return arms


def family_constraints(
    games: Sequence[Game], edges: Sequence[Edge], n: int
) -> tuple[
    list[list[F]],
    list[F],
    list[list[F]],
    list[F],
    list[tuple[Any, ...]],
]:
    grand = (1 << n) - 1
    variable_count = len(games) * n
    inequalities: list[list[F]] = []
    inequality_rhs: list[F] = []
    equalities: list[list[F]] = []
    equality_rhs: list[F] = []
    metadata: list[tuple[Any, ...]] = []

    for game_index, game in enumerate(games):
        row = [F(0)] * variable_count
        for player in range(n):
            row[game_index * n + player] = 1
        equalities.append(row)
        equality_rhs.append(game[grand])

        for coalition in range(1, grand):
            row = [F(0)] * variable_count
            for player in members(coalition, n):
                row[game_index * n + player] = -1
            inequalities.append(row)
            inequality_rhs.append(-game[coalition])
            metadata.append(("core", game_index, coalition))

    for lower, upper, coalition in edges:
        for player in members(coalition, n):
            row = [F(0)] * variable_count
            row[lower * n + player] = 1
            row[upper * n + player] = -1
            inequalities.append(row)
            inequality_rhs.append(F(0))
            metadata.append(
                ("monotonicity", lower, upper, coalition, player)
            )
    return (
        inequalities,
        inequality_rhs,
        equalities,
        equality_rhs,
        metadata,
    )


def box_family_constraints(
    games: Sequence[Game], edges: Sequence[Edge], n: int
) -> tuple[
    list[list[F]],
    list[F],
    list[list[F]],
    list[F],
    list[tuple[Any, ...]],
]:
    """Family constraints using only exact coordinate projections."""

    grand = (1 << n) - 1
    variable_count = len(games) * n
    inequalities: list[list[F]] = []
    inequality_rhs: list[F] = []
    equalities: list[list[F]] = []
    equality_rhs: list[F] = []
    metadata: list[tuple[Any, ...]] = []

    for game_index, game in enumerate(games):
        row = [F(0)] * variable_count
        for player in range(n):
            row[game_index * n + player] = F(1)
        equalities.append(row)
        equality_rhs.append(game[grand])

        for player in range(n):
            lower = [F(0)] * variable_count
            lower[game_index * n + player] = F(-1)
            inequalities.append(lower)
            inequality_rhs.append(-game[1 << player])
            metadata.append(("box_lower", game_index, player))

            upper = [F(0)] * variable_count
            upper[game_index * n + player] = F(1)
            inequalities.append(upper)
            inequality_rhs.append(
                game[grand] - game[grand ^ (1 << player)]
            )
            metadata.append(("box_upper", game_index, player))

    for lower, upper, coalition in edges:
        for player in members(coalition, n):
            row = [F(0)] * variable_count
            row[lower * n + player] = F(1)
            row[upper * n + player] = F(-1)
            inequalities.append(row)
            inequality_rhs.append(F(0))
            metadata.append(
                ("monotonicity", lower, upper, coalition, player)
            )
    return (
        inequalities,
        inequality_rhs,
        equalities,
        equality_rhs,
        metadata,
    )


def _family_slack_float_from_constraints(
    games: Sequence[Game],
    n: int,
    constraints: tuple[
        list[list[F]],
        list[F],
        list[list[F]],
        list[F],
        list[tuple[Any, ...]],
    ],
) -> tuple[float | None, list[list[float]] | None]:
    (
        inequalities,
        inequality_rhs,
        equalities,
        equality_rhs,
        metadata,
    ) = constraints
    allocation_variables = len(games) * n
    variable_count = allocation_variables + 1
    slack_index = allocation_variables

    float_inequalities: list[list[float]] = []
    float_rhs = [float(value) for value in inequality_rhs]
    for row, tag in zip(inequalities, metadata, strict=True):
        extended = [float(value) for value in row] + [0.0]
        if tag[0] == "monotonicity":
            extended[slack_index] = 1.0
        float_inequalities.append(extended)
    float_equalities = [
        [float(value) for value in row] + [0.0]
        for row in equalities
    ]
    objective = [0.0] * variable_count
    objective[slack_index] = -1.0
    result = scipy_linprog(
        objective,
        A_ub=float_inequalities,
        b_ub=float_rhs,
        A_eq=float_equalities,
        b_eq=[float(value) for value in equality_rhs],
        bounds=[(None, None)] * variable_count,
        method="highs",
    )
    if not result.success:
        return None, None
    allocations = [
        [
            float(result.x[game_index * n + player])
            for player in range(n)
        ]
        for game_index in range(len(games))
    ]
    return float(result.x[slack_index]), allocations


def family_slack_float(
    games: Sequence[Game], edges: Sequence[Edge], n: int
) -> tuple[float | None, list[list[float]] | None]:
    """Maximize the common monotonicity margin.

    A nonnegative optimum is equivalent to family feasibility.  A negative
    optimum proves float-screened infeasibility of the weak monotonicity
    constraints.
    """

    return _family_slack_float_from_constraints(
        games, n, family_constraints(games, edges, n)
    )


def box_family_slack_float(
    games: Sequence[Game], edges: Sequence[Edge], n: int
) -> tuple[float | None, list[list[float]] | None]:
    """Maximize the margin over coordinate-box sections of the cores."""

    return _family_slack_float_from_constraints(
        games, n, box_family_constraints(games, edges, n)
    )


def _family_slack_exact_from_constraints(
    games: Sequence[Game],
    n: int,
    constraints: tuple[
        list[list[F]],
        list[F],
        list[list[F]],
        list[F],
        list[tuple[Any, ...]],
    ],
) -> tuple[F | None, list[Point] | None]:
    inequalities, rhs, equalities, eq_rhs, metadata = constraints
    allocation_variables = len(games) * n
    extended_inequalities: list[list[F]] = []
    for row, tag in zip(inequalities, metadata, strict=True):
        extended = list(row) + [F(0)]
        if tag[0] == "monotonicity":
            extended[-1] = 1
        extended_inequalities.append(extended)
    extended_equalities = [list(row) + [F(0)] for row in equalities]
    objective = [F(0)] * (allocation_variables + 1)
    objective[-1] = -1
    status, optimum, point = exact_linprog(
        objective,
        extended_inequalities,
        rhs,
        extended_equalities,
        eq_rhs,
    )
    if status != "optimal" or optimum is None or point is None:
        return None, None
    allocations = [
        tuple(point[index * n : (index + 1) * n])
        for index in range(len(games))
    ]
    return -optimum, allocations


def family_slack_exact(
    games: Sequence[Game], edges: Sequence[Edge], n: int
) -> tuple[F | None, list[Point] | None]:
    return _family_slack_exact_from_constraints(
        games, n, family_constraints(games, edges, n)
    )


def box_family_slack_exact(
    games: Sequence[Game], edges: Sequence[Edge], n: int
) -> tuple[F | None, list[Point] | None]:
    return _family_slack_exact_from_constraints(
        games, n, box_family_constraints(games, edges, n)
    )


def common_core_gap_float(
    games: Sequence[Game], n: int
) -> float | None:
    """Return the minimum common-core budget minus the games' common total.

    A positive value means that the cores have empty intersection.  When the
    value is nonpositive, one common allocation can be selected for every game,
    so the family is automatically compatible with coalitional monotonicity.
    """

    grand = (1 << n) - 1
    assert all(game[grand] == games[0][grand] for game in games)
    envelope = [
        max(game[coalition] for game in games)
        for coalition in range(1 << n)
    ]
    rows: list[list[float]] = []
    rhs: list[float] = []
    for coalition in range(1, grand):
        rows.append(
            [
                -1.0 if coalition >> player & 1 else 0.0
                for player in range(n)
            ]
        )
        rhs.append(-float(envelope[coalition]))
    result = scipy_linprog(
        [1.0] * n,
        A_ub=rows,
        b_ub=rhs,
        bounds=[(None, None)] * n,
        method="highs",
    )
    if not result.success:
        return None
    return float(result.fun) - float(games[0][grand])


def common_core_gap_exact(games: Sequence[Game], n: int) -> F | None:
    grand = (1 << n) - 1
    assert all(game[grand] == games[0][grand] for game in games)
    envelope = [
        max(game[coalition] for game in games)
        for coalition in range(1 << n)
    ]
    inequalities: list[list[F]] = []
    rhs: list[F] = []
    for coalition in range(1, grand):
        inequalities.append(
            [
                F(-1) if coalition >> player & 1 else F(0)
                for player in range(n)
            ]
        )
        rhs.append(-envelope[coalition])
    status, optimum, _ = exact_linprog(
        [F(1)] * n, inequalities, rhs
    )
    if status != "optimal" or optimum is None:
        return None
    return optimum - games[0][grand]


def minimum_coordinate_width(game: Game, n: int) -> F:
    """Minimum range of a player's payoff over an exact game's core."""

    grand = (1 << n) - 1
    return min(
        game[grand]
        - game[grand ^ (1 << player)]
        - game[1 << player]
        for player in range(n)
    )


def protected_path_width_bound(
    family: "Family", n: int
) -> tuple[F, int, list[int]] | None:
    """Best coordinate-width upper bound on a uniform edge margin."""

    if not family.edges:
        return None
    grand = (1 << n) - 1
    order = sorted(
        range(len(family.games)),
        key=lambda node: (family.depths[node], node),
    )
    best: tuple[F, int, list[int]] | None = None
    for player in range(n):
        adjacency: dict[int, list[int]] = {}
        for lower, upper, coalition in family.edges:
            if coalition >> player & 1:
                adjacency.setdefault(lower, []).append(upper)
        for start in order:
            distance = {start: 0}
            paths = {start: [start]}
            for node in order:
                if node not in distance:
                    continue
                for child in adjacency.get(node, []):
                    candidate = distance[node] + 1
                    if candidate > distance.get(child, -1):
                        distance[child] = candidate
                        paths[child] = paths[node] + [child]
            lower_bound = family.games[start][1 << player]
            for end, length in distance.items():
                if length == 0:
                    continue
                upper_bound = (
                    family.games[end][grand]
                    - family.games[end][grand ^ (1 << player)]
                )
                bound = (upper_bound - lower_bound) / length
                candidate = (bound, player, paths[end])
                if best is None or candidate[0] < best[0]:
                    best = candidate
    return best


def family_feasible_exact(
    games: Sequence[Game], edges: Sequence[Edge], n: int
) -> tuple[bool, list[Point] | None]:
    inequalities, rhs, equalities, eq_rhs, _ = family_constraints(
        games, edges, n
    )
    variable_count = len(games) * n
    status, _, point = exact_linprog(
        [F(0)] * variable_count,
        inequalities,
        rhs,
        equalities,
        eq_rhs,
    )
    if status != "optimal" or point is None:
        return False, None
    allocations = [
        tuple(point[game_index * n : (game_index + 1) * n])
        for game_index in range(len(games))
    ]
    return True, allocations


def verify_family_allocations(
    games: Sequence[Game],
    edges: Sequence[Edge],
    allocations: Sequence[Point],
    n: int,
    margin: F = F(0),
) -> bool:
    grand = (1 << n) - 1
    if len(allocations) != len(games):
        return False
    for game, point in zip(games, allocations, strict=True):
        if sum(point) != game[grand]:
            return False
        if any(
            coalition_sum(point, coalition) < game[coalition]
            for coalition in range(1, grand)
        ):
            return False
    for lower, upper, coalition in edges:
        if any(
            allocations[upper][player]
            < allocations[lower][player] + margin
            for player in members(coalition, n)
        ):
            return False
    return True


def reconstructed_family_witness(
    games: Sequence[Game],
    edges: Sequence[Edge],
    n: int,
    max_denominator: int = 1_000_000,
) -> tuple[F, list[Point]] | None:
    slack, float_allocations = family_slack_float(games, edges, n)
    if slack is None or float_allocations is None:
        return None
    margin = F(slack).limit_denominator(max_denominator)
    allocations = [
        tuple(
            F(value).limit_denominator(max_denominator)
            for value in point
        )
        for point in float_allocations
    ]
    if not verify_family_allocations(
        games, edges, allocations, n, margin
    ):
        return None
    return margin, allocations


def _reconstructed_family_optimum_from_constraints(
    games: Sequence[Game],
    n: int,
    constraints: tuple[
        list[list[F]],
        list[F],
        list[list[F]],
        list[F],
        list[tuple[Any, ...]],
    ],
    max_denominator: int = 1_000_000,
) -> tuple[F, list[Point], dict[str, Any]] | None:
    """Reconstruct and exactly verify primal and dual margin optima."""

    inequalities, rhs, equalities, eq_rhs, metadata = constraints
    allocation_variables = len(games) * n
    variable_count = allocation_variables + 1
    extended_inequalities: list[list[F]] = []
    for row, tag in zip(inequalities, metadata, strict=True):
        extended = list(row) + [F(0)]
        if tag[0] == "monotonicity":
            extended[-1] = F(1)
        extended_inequalities.append(extended)
    extended_equalities = [list(row) + [F(0)] for row in equalities]
    objective = [F(0)] * variable_count
    objective[-1] = F(-1)

    result = scipy_linprog(
        [float(value) for value in objective],
        A_ub=[
            [float(value) for value in row]
            for row in extended_inequalities
        ],
        b_ub=[float(value) for value in rhs],
        A_eq=[
            [float(value) for value in row]
            for row in extended_equalities
        ],
        b_eq=[float(value) for value in eq_rhs],
        bounds=[(None, None)] * variable_count,
        method="highs",
    )
    if not result.success:
        return None

    primal = [
        F(value).limit_denominator(max_denominator)
        for value in result.x
    ]
    margin = primal[-1]
    allocations = [
        tuple(primal[index * n : (index + 1) * n])
        for index in range(len(games))
    ]
    if any(
        exact_dot(row, primal) > bound
        for row, bound in zip(
            extended_inequalities, rhs, strict=True
        )
    ):
        return None
    if any(
        exact_dot(row, primal) != bound
        for row, bound in zip(
            extended_equalities, eq_rhs, strict=True
        )
    ):
        return None

    inequality_weights = [
        F(value).limit_denominator(max_denominator)
        for value in result.ineqlin.marginals
    ]
    equality_weights = [
        F(value).limit_denominator(max_denominator)
        for value in result.eqlin.marginals
    ]
    if any(weight > 0 for weight in inequality_weights):
        return None
    active_inequality_rows = [
        (weight, row)
        for weight, row in zip(
            inequality_weights,
            extended_inequalities,
            strict=True,
        )
        if weight
    ]
    active_equality_rows = [
        (weight, row)
        for weight, row in zip(
            equality_weights,
            extended_equalities,
            strict=True,
        )
        if weight
    ]
    for column in range(variable_count):
        dual_column = sum(
            (
                weight * row[column]
                for weight, row in active_inequality_rows
            ),
            start=F(0),
        )
        dual_column += sum(
            (
                weight * row[column]
                for weight, row in active_equality_rows
            ),
            start=F(0),
        )
        if dual_column != objective[column]:
            return None
    dual_objective = sum(
        (
            weight * bound
            for weight, bound in zip(
                inequality_weights, rhs, strict=True
            )
        ),
        start=F(0),
    )
    dual_objective += sum(
        (
            weight * bound
            for weight, bound in zip(
                equality_weights, eq_rhs, strict=True
            )
        ),
        start=F(0),
    )
    if dual_objective != -margin:
        return None

    certificate = {
        "objective": str(dual_objective),
        "active_inequalities": [
            {
                "row": metadata[index],
                "weight": str(weight),
            }
            for index, weight in enumerate(inequality_weights)
            if weight
        ],
        "active_equalities": [
            {"game": index, "weight": str(weight)}
            for index, weight in enumerate(equality_weights)
            if weight
        ],
    }
    return margin, allocations, certificate


def _active_basis_family_optimum_from_constraints(
    games: Sequence[Game],
    n: int,
    constraints: tuple[
        list[list[F]],
        list[F],
        list[list[F]],
        list[F],
        list[tuple[Any, ...]],
    ],
) -> tuple[F, list[Point], dict[str, Any]] | None:
    """Reconstruct a degenerate rational optimum from a HiGHS active basis."""

    import numpy as np
    from scipy.linalg import qr
    from sympy import Rational
    from sympy.polys.matrices import DomainMatrix

    inequalities, rhs, equalities, eq_rhs, metadata = constraints
    allocation_variables = len(games) * n
    variable_count = allocation_variables + 1
    extended_inequalities = []
    for row, tag in zip(inequalities, metadata, strict=True):
        extended = list(row) + [F(0)]
        if tag[0] == "monotonicity":
            extended[-1] = F(1)
        extended_inequalities.append(extended)
    extended_equalities = [list(row) + [F(0)] for row in equalities]
    float_inequalities = np.asarray(
        [
            [float(value) for value in row]
            for row in extended_inequalities
        ]
    )
    float_equalities = np.asarray(
        [
            [float(value) for value in row]
            for row in extended_equalities
        ]
    )
    objective = np.zeros(variable_count)
    objective[-1] = -1
    result = scipy_linprog(
        objective,
        A_ub=float_inequalities,
        b_ub=np.asarray([float(value) for value in rhs]),
        A_eq=float_equalities,
        b_eq=np.asarray([float(value) for value in eq_rhs]),
        bounds=[(None, None)] * variable_count,
        method="highs",
    )
    if not result.success:
        return None

    active = np.flatnonzero(result.ineqlin.residual < 1e-8)
    support = np.flatnonzero(
        abs(result.ineqlin.marginals) > 1e-10
    )
    support_set = set(int(index) for index in support)
    base_float = np.vstack(
        [float_equalities, float_inequalities[support]]
    )
    base_rank = int(np.linalg.matrix_rank(base_float))
    if base_rank != len(base_float):
        return None
    q_matrix, _ = qr(base_float.T, mode="full")
    null_basis = q_matrix[:, base_rank:]
    remaining = np.asarray(
        [
            int(index)
            for index in active
            if int(index) not in support_set
        ]
    )
    needed = variable_count - base_rank
    if needed < 0 or len(remaining) < needed:
        return None
    projected = float_inequalities[remaining] @ null_basis
    _, _, pivots = qr(
        projected.T, pivoting=True, mode="economic"
    )
    completion = remaining[pivots[:needed]]
    basis_rows = (
        extended_equalities
        + [
            extended_inequalities[int(index)]
            for index in support
        ]
        + [
            extended_inequalities[int(index)]
            for index in completion
        ]
    )
    basis_rhs = (
        eq_rhs
        + [rhs[int(index)] for index in support]
        + [rhs[int(index)] for index in completion]
    )
    if len(basis_rows) != variable_count:
        return None
    if (
        np.linalg.matrix_rank(
            np.asarray(
                [
                    [float(value) for value in row]
                    for row in basis_rows
                ]
            )
        )
        != variable_count
    ):
        return None

    def domain_matrix(
        rows: Sequence[Sequence[F]],
    ) -> DomainMatrix:
        values: dict[int, dict[int, Any]] = {}
        for row_index, row in enumerate(rows):
            for column_index, value in enumerate(row):
                if value:
                    values.setdefault(row_index, {})[
                        column_index
                    ] = Rational(
                        value.numerator, value.denominator
                    )
        return DomainMatrix.from_dict_sympy(
            len(rows), len(rows[0]), values
        ).to_field()

    try:
        basis = domain_matrix(basis_rows)
        solution_numerator, solution_denominator = basis.solve_den(
            domain_matrix([[value] for value in basis_rhs]),
            method="rref",
        )
        dual_rhs = [[F(0)] for _ in range(variable_count)]
        dual_rhs[-1][0] = F(-1)
        dual_numerator, dual_denominator = (
            basis.transpose().solve_den(
                domain_matrix(dual_rhs), method="rref"
            )
        )
    except Exception:
        return None

    solution_matrix = solution_numerator.to_Matrix()
    solution_scale = F(str(solution_denominator))
    primal = [
        F(str(solution_matrix[index, 0])) / solution_scale
        for index in range(variable_count)
    ]
    if any(
        exact_dot(row, primal) > bound
        for row, bound in zip(
            extended_inequalities, rhs, strict=True
        )
    ):
        return None
    if any(
        exact_dot(row, primal) != bound
        for row, bound in zip(
            extended_equalities, eq_rhs, strict=True
        )
    ):
        return None

    dual_matrix = dual_numerator.to_Matrix()
    dual_scale = F(str(dual_denominator))
    basis_dual = [
        F(str(dual_matrix[index, 0])) / dual_scale
        for index in range(variable_count)
    ]
    equality_weights = basis_dual[: len(equalities)]
    inequality_weights = [F(0)] * len(inequalities)
    for index, weight in zip(
        support,
        basis_dual[
            len(equalities) : len(equalities) + len(support)
        ],
        strict=True,
    ):
        inequality_weights[int(index)] = weight
    for index, weight in zip(
        completion,
        basis_dual[len(equalities) + len(support) :],
        strict=True,
    ):
        inequality_weights[int(index)] = weight
    if any(weight > 0 for weight in inequality_weights):
        return None
    margin = primal[-1]
    dual_objective = sum(
        (
            weight * bound
            for weight, bound in zip(
                inequality_weights, rhs, strict=True
            )
        ),
        start=F(0),
    )
    dual_objective += sum(
        (
            weight * bound
            for weight, bound in zip(
                equality_weights, eq_rhs, strict=True
            )
        ),
        start=F(0),
    )
    if dual_objective != -margin:
        return None

    allocations = [
        tuple(
            primal[index * n : (index + 1) * n]
        )
        for index in range(len(games))
    ]
    certificate = {
        "objective": str(dual_objective),
        "active_inequalities": [
            {
                "row": metadata[index],
                "weight": str(weight),
            }
            for index, weight in enumerate(inequality_weights)
            if weight
        ],
        "active_equalities": [
            {"game": index, "weight": str(weight)}
            for index, weight in enumerate(equality_weights)
            if weight
        ],
        "reconstruction": "exact_active_basis",
    }
    return margin, allocations, certificate


def reconstructed_family_optimum(
    games: Sequence[Game],
    edges: Sequence[Edge],
    n: int,
    max_denominator: int = 1_000_000,
) -> tuple[F, list[Point], dict[str, Any]] | None:
    constraints = family_constraints(games, edges, n)
    optimum = _reconstructed_family_optimum_from_constraints(
        games,
        n,
        constraints,
        max_denominator,
    )
    if optimum is not None:
        return optimum
    return _active_basis_family_optimum_from_constraints(
        games, n, constraints
    )


def reconstructed_box_family_optimum(
    games: Sequence[Game],
    edges: Sequence[Edge],
    n: int,
    max_denominator: int = 1_000_000,
) -> tuple[F, list[Point], dict[str, Any]] | None:
    constraints = box_family_constraints(games, edges, n)
    optimum = _reconstructed_family_optimum_from_constraints(
        games,
        n,
        constraints,
        max_denominator,
    )
    if optimum is not None:
        return optimum
    return _active_basis_family_optimum_from_constraints(
        games, n, constraints
    )


def farkas_certificate(
    games: Sequence[Game], edges: Sequence[Edge], n: int
) -> dict[str, Any] | None:
    inequalities, rhs, equalities, eq_rhs, metadata = family_constraints(
        games, edges, n
    )
    inequality_count = len(inequalities)
    equality_count = len(equalities)
    allocation_variables = len(games) * n
    certificate_variables = inequality_count + equality_count

    certificate_equalities: list[list[F]] = []
    certificate_rhs: list[F] = []
    for column in range(allocation_variables):
        certificate_equalities.append(
            [row[column] for row in inequalities]
            + [row[column] for row in equalities]
        )
        certificate_rhs.append(F(0))
    certificate_equalities.append(list(rhs) + list(eq_rhs))
    certificate_rhs.append(F(-1))

    nonnegativity: list[list[F]] = []
    nonnegativity_rhs: list[F] = []
    for index in range(inequality_count):
        row = [F(0)] * certificate_variables
        row[index] = -1
        nonnegativity.append(row)
        nonnegativity_rhs.append(F(0))

    status, _, point = exact_linprog(
        [F(0)] * certificate_variables,
        nonnegativity,
        nonnegativity_rhs,
        certificate_equalities,
        certificate_rhs,
    )
    if status != "optimal" or point is None:
        return None
    inequality_weights = point[:inequality_count]
    equality_weights = point[inequality_count:]
    if any(weight < 0 for weight in inequality_weights):
        return None
    for column in range(allocation_variables):
        total = sum(
            (
                inequality_weights[index]
                * inequalities[index][column]
                for index in range(inequality_count)
            ),
            start=F(0),
        )
        total += sum(
            (
                equality_weights[index] * equalities[index][column]
                for index in range(equality_count)
            ),
            start=F(0),
        )
        if total != 0:
            return None
    contradiction = sum(
        (
            inequality_weights[index] * rhs[index]
            for index in range(inequality_count)
        ),
        start=F(0),
    )
    contradiction += sum(
        (
            equality_weights[index] * eq_rhs[index]
            for index in range(equality_count)
        ),
        start=F(0),
    )
    if contradiction >= 0:
        return None
    return {
        "rhs": str(contradiction),
        "active_inequalities": [
            {
                "row": metadata[index],
                "weight": str(weight),
            }
            for index, weight in enumerate(inequality_weights)
            if weight
        ],
        "active_equalities": [
            {"game": index, "weight": str(weight)}
            for index, weight in enumerate(equality_weights)
            if weight
        ],
    }


@dataclass
class Family:
    games: list[Game]
    edges: list[Edge]
    depths: list[int]
    shape: str


def star_family(
    base: Game,
    arms: Sequence[tuple[int, F, Game]],
    n: int,
    max_nodes: int,
    shape: str,
) -> Family:
    games = [base]
    edges: list[Edge] = []
    for coalition, _, game in arms[: max_nodes - 1]:
        edges.append((0, len(games), coalition))
        games.append(game)
    return Family(games, edges, [0] + [1] * len(edges), shape)


def grow_family(
    base: Game,
    n: int,
    args: argparse.Namespace,
    rng: random.Random,
    base_arms: Sequence[tuple[int, F, Game]] | None = None,
) -> Family:
    games = [base]
    edges: list[Edge] = []
    depths = [0]
    game_index = {base: 0}
    frontier = [0]

    for depth in range(args.depth):
        next_frontier: list[int] = []
        rng.shuffle(frontier)
        for parent in frontier:
            if depth == 0 and parent == 0 and base_arms is not None:
                arms = list(base_arms)
                rng.shuffle(arms)
            else:
                arms = arm_pool(
                    games[parent],
                    n,
                    F(1, args.bump_denominator),
                    args.bump_steps,
                    args.coalitions_per_node,
                    args.require_monotone,
                    args.require_superadditive,
                    rng,
                )
            if args.dag_fraction_denominator > 1:
                arms = [
                    (
                        coalition,
                        delta / args.dag_fraction_denominator,
                        bump_game(
                            games[parent],
                            coalition,
                            delta / args.dag_fraction_denominator,
                        ),
                    )
                    for coalition, delta, _ in arms
                ]
            elif args.dag_backoff_steps:
                backoff = (
                    F(args.dag_backoff_steps, args.bump_denominator)
                )
                arms = [
                    (
                        coalition,
                        (
                            delta - backoff
                            if delta > backoff
                            else delta / 2
                        ),
                        bump_game(
                            games[parent],
                            coalition,
                            (
                                delta - backoff
                                if delta > backoff
                                else delta / 2
                            ),
                        ),
                    )
                    for coalition, delta, _ in arms
                ]
            if args.require_positive_width:
                arms = [
                    arm
                    for arm in arms
                    if minimum_coordinate_width(arm[2], n) > 0
                ]
            added = 0
            for coalition, _, child_game in arms:
                if child_game in game_index:
                    child = game_index[child_game]
                elif len(games) < args.max_nodes:
                    child = len(games)
                    game_index[child_game] = child
                    games.append(child_game)
                    depths.append(depth + 1)
                    next_frontier.append(child)
                else:
                    continue
                edge = (parent, child, coalition)
                if edge not in edges:
                    edges.append(edge)
                added += 1
                if added >= args.branching:
                    break
        frontier = next_frontier
        if not frontier or len(games) >= args.max_nodes:
            break
    return Family(games, edges, depths, "bump_dag")


def serial_game(game: Game) -> list[str]:
    return [str(value) for value in game]


def serial_family(
    family: Family,
    slack: float | None,
    common_core_gap: float | None = None,
    exact_verified: bool | None = None,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "shape": family.shape,
        "games": [serial_game(game) for game in family.games],
        "edges": [
            {
                "lower": lower,
                "upper": upper,
                "coalition": coalition,
                "members": [
                    player + 1
                    for player in members(
                        coalition,
                        (len(family.games[0]) - 1).bit_length(),
                    )
                ],
                "delta": str(
                    family.games[upper][coalition]
                    - family.games[lower][coalition]
                ),
            }
            for lower, upper, coalition in family.edges
        ],
        "depths": family.depths,
        "max_min_monotonicity_margin_float": slack,
    }
    if common_core_gap is not None:
        payload["common_core_budget_gap_float"] = common_core_gap
    n = (len(family.games[0]) - 1).bit_length()
    minimum_width = min(
        minimum_coordinate_width(game, n) for game in family.games
    )
    payload["minimum_coordinate_width"] = str(minimum_width)
    path_bound = protected_path_width_bound(family, n)
    if path_bound is not None:
        bound, player, path = path_bound
        payload["protected_path_width_bound"] = {
            "bound": str(bound),
            "player": player + 1,
            "path": path,
        }
    if slack is not None and minimum_width > 0:
        payload["margin_to_minimum_width_ratio_float"] = (
            slack / float(minimum_width)
        )
    if exact_verified is not None:
        payload["all_games_exact_verified"] = exact_verified
    return payload


def exact_verify_family(family: Family, n: int) -> bool:
    grand = (1 << n) - 1
    for lower, upper, coalition in family.edges:
        if (
            lower not in range(len(family.games))
            or upper not in range(len(family.games))
            or coalition not in range(1, grand + 1)
        ):
            return False
        lower_game = family.games[lower]
        upper_game = family.games[upper]
        if upper_game[coalition] <= lower_game[coalition]:
            return False
        if any(
            upper_game[target] != lower_game[target]
            for target in range(grand + 1)
            if target != coalition
        ):
            return False
    for game in family.games:
        certificate = reconstructed_exactness_certificate(game, n)
        if certificate is None:
            certificate = exactness_certificate(game, n)
        if certificate is None:
            return False
    return True


def search(args: argparse.Namespace) -> int:
    rng = random.Random(args.seed)
    started = time.time()
    best_family: Family | None = None
    best_slack: float | None = None
    best_score: float | None = None
    best_common_core_gap: float | None = None
    families = 0
    games_screened = 0
    empty_common_core_families = 0
    exact_hits = 0
    shape_counts: dict[str, int] = {}
    nontrivial_shape_counts: dict[str, int] = {}

    for trial in range(1, args.trials + 1):
        total = rng.randint(args.n, args.total_max)
        style_draw = rng.random()
        if style_draw < args.simplex_probability:
            points = simplex_points(args.n, total)
        elif style_draw < (
            args.simplex_probability + args.near_simplex_probability
        ):
            points = near_simplex_points(args.n, total, rng)
        else:
            point_count = rng.randint(args.min_points, args.max_points)
            points = [
                random_composition(args.n, total, rng)
                for _ in range(point_count)
            ]
        base = game_from_points(points)
        base_arms = arm_pool(
            base,
            args.n,
            F(1, args.bump_denominator),
            args.bump_steps,
            0,
            args.require_monotone,
            args.require_superadditive,
            rng,
        )
        if not base_arms:
            continue
        cosingletons = [
            arm
            for arm in base_arms
            if len(members(arm[0], args.n)) == args.n - 1
        ]
        dense_order = list(cosingletons)
        dense_order.extend(
            arm for arm in base_arms if arm not in cosingletons
        )
        if args.require_positive_width:
            dense_order = [
                arm
                for arm in dense_order
                if minimum_coordinate_width(arm[2], args.n) > 0
            ]
            cosingletons = [
                arm
                for arm in cosingletons
                if minimum_coordinate_width(arm[2], args.n) > 0
            ]
        candidates = [
            grow_family(base, args.n, args, rng, base_arms),
            star_family(
                base,
                dense_order,
                args.n,
                args.max_nodes,
                "dense_star",
            ),
        ]
        if len(cosingletons) >= 2:
            candidates.append(
                star_family(
                    base,
                    cosingletons,
                    args.n,
                    args.max_nodes,
                    "cosingleton_star",
                )
            )

        for family in candidates:
            if not family.edges:
                continue
            families += 1
            shape_counts[family.shape] = (
                shape_counts.get(family.shape, 0) + 1
            )
            games_screened += len(family.games)
            common_core_gap = common_core_gap_float(
                family.games, args.n
            )
            if common_core_gap is None or common_core_gap <= 1e-8:
                continue
            empty_common_core_families += 1
            nontrivial_shape_counts[family.shape] = (
                nontrivial_shape_counts.get(family.shape, 0) + 1
            )
            slack, _ = family_slack_float(
                family.games, family.edges, args.n
            )
            if slack is None:
                continue
            minimum_width = min(
                minimum_coordinate_width(game, args.n)
                for game in family.games
            )
            score = (
                slack / float(minimum_width)
                if minimum_width > 0
                else slack
            )
            if best_score is None or score < best_score:
                best_score = score
                best_slack = slack
                best_family = family
                best_common_core_gap = common_core_gap
                print(
                    f"best trial={trial} shape={family.shape} "
                    f"nodes={len(family.games)} "
                    f"edges={len(family.edges)} margin={slack:.9g} "
                    f"score={score:.9g} common_gap={common_core_gap}",
                    flush=True,
                )

            if slack < -args.infeasible_tolerance:
                all_exact = exact_verify_family(family, args.n)
                if not all_exact:
                    continue
                feasible, _ = family_feasible_exact(
                    family.games, family.edges, args.n
                )
                if feasible:
                    continue
                certificate = farkas_certificate(
                    family.games, family.edges, args.n
                )
                if certificate is None:
                    continue
                exact_hits += 1
                payload = {
                    "status": "certified_infeasible",
                    "claim": (
                        "No coalitionally monotonic core selection exists "
                        "on any domain containing this finite exact-game "
                        "family."
                    ),
                    "n": args.n,
                    "seed": args.seed,
                    "trial": trial,
                    "family": serial_family(
                        family, slack, common_core_gap, True
                    ),
                    "farkas_certificate": certificate,
                    "elapsed_seconds": time.time() - started,
                }
                Path(args.output).write_text(
                    json.dumps(payload, indent=2) + "\n"
                )
                print(json.dumps(payload, indent=2), flush=True)
                return 1

        if trial % args.progress_every == 0:
            print(
                f"trial={trial} families={families} "
                f"nontrivial={empty_common_core_families} "
                f"games={games_screened} best_score={best_score} "
                f"elapsed={time.time() - started:.1f}s",
                flush=True,
            )

    best_exact = (
        exact_verify_family(best_family, args.n)
        if best_family is not None and args.verify_best
        else None
    )
    summary = {
        "status": "no_infeasible_family_found",
        "n": args.n,
        "seed": args.seed,
        "trials": args.trials,
        "families_screened": families,
        "games_screened": games_screened,
        "empty_common_core_families": empty_common_core_families,
        "shape_counts": shape_counts,
        "nontrivial_shape_counts": nontrivial_shape_counts,
        "certified_hits": exact_hits,
        "best_normalized_score": best_score,
        "best_family": (
            serial_family(
                best_family,
                best_slack,
                best_common_core_gap,
                best_exact,
            )
            if best_family is not None
            else None
        ),
        "search_configuration": {
            key: value
            for key, value in vars(args).items()
            if key != "output"
        },
        "elapsed_seconds": time.time() - started,
    }
    Path(args.output).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--min-points", type=int, default=2)
    parser.add_argument("--max-points", type=int, default=6)
    parser.add_argument("--total-max", type=int, default=16)
    parser.add_argument("--simplex-probability", type=float, default=0.2)
    parser.add_argument(
        "--near-simplex-probability", type=float, default=0.3
    )
    parser.add_argument("--depth", type=int, default=2)
    parser.add_argument("--branching", type=int, default=3)
    parser.add_argument("--max-nodes", type=int, default=16)
    parser.add_argument("--coalitions-per-node", type=int, default=30)
    parser.add_argument("--bump-denominator", type=int, default=2)
    parser.add_argument("--bump-steps", type=int, default=6)
    parser.add_argument("--dag-backoff-steps", type=int, default=1)
    parser.add_argument("--dag-fraction-denominator", type=int, default=1)
    parser.add_argument("--require-monotone", action="store_true")
    parser.add_argument("--require-superadditive", action="store_true")
    parser.add_argument("--require-positive-width", action="store_true")
    parser.add_argument("--verify-best", action="store_true")
    parser.add_argument("--infeasible-tolerance", type=float, default=1e-7)
    parser.add_argument("--progress-every", type=int, default=10)
    parser.add_argument(
        "--output", default="results/family_search.json"
    )
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(search(parse_args()))
