#!/usr/bin/env python3
"""Exactly verify the projected-Shapley monotonicity counterexample."""

from __future__ import annotations

import math
import sys
from fractions import Fraction as F


sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_facet_search import load_facets


N = 5
GRAND = 31
CHANGED = 7
LOWER = tuple(
    F(value)
    for value in (
        "0 0 0 333/1000 0 333/1000 0 333/1000 "
        "0 333/1000 0 333/1000 0 333/1000 0 167/500 "
        "0 0 0 333/1000 0 167/500 0 167/500 "
        "333/1000 333/1000 333/500 333/500 333/500 "
        "667/1000 333/500 1"
    ).split()
)
UPPER = tuple(
    value + (F(1, 1000) if coalition == CHANGED else 0)
    for coalition, value in enumerate(LOWER)
)
LOWER_PROJECTION = (
    F(8819, 42000),
    F(5167, 42000),
    F(5167, 42000),
    F(8003, 28000),
    F(4337, 16800),
)
UPPER_PROJECTION = (
    F(4409, 21000),
    F(323, 2625),
    F(323, 2625),
    F(3001, 10500),
    F(1807, 7000),
)
ACTIVE = (3, 5)
MULTIPLIERS = (
    (F(311, 84000), F(283, 84000)),
    (F(307, 84000), F(93, 28000)),
)


def shapley(game: tuple[F, ...]) -> tuple[F, ...]:
    values = []
    for player in range(N):
        value = F(0)
        for coalition in range(GRAND + 1):
            if coalition >> player & 1:
                continue
            size = coalition.bit_count()
            weight = F(
                math.factorial(size) * math.factorial(N - size - 1),
                math.factorial(N),
            )
            value += weight * (
                game[coalition | (1 << player)] - game[coalition]
            )
        values.append(value)
    return tuple(values)


def normal(coalition: int) -> tuple[F, ...]:
    size = coalition.bit_count()
    return tuple(
        F(N * int(bool(coalition >> player & 1)) - size, N)
        for player in range(N)
    )


def verify_game(game: tuple[F, ...]) -> None:
    facets, _ = load_facets()
    assert game[0] == 0 and game[GRAND] == 1
    for coalition in range(GRAND + 1):
        for player in range(N):
            if not coalition >> player & 1:
                assert game[coalition] <= game[coalition | (1 << player)]
    assert all(
        sum(F(coefficient) * game[coalition] for coalition, coefficient in enumerate(facet))
        >= 0
        for facet in facets
    )


def verify_projection(
    game: tuple[F, ...],
    projection: tuple[F, ...],
    multipliers: tuple[F, ...],
) -> None:
    reference = shapley(game)
    assert sum(projection) == game[GRAND]
    assert all(
        sum(projection[player] for player in range(N) if coalition >> player & 1)
        >= game[coalition]
        for coalition in range(1, GRAND)
    )
    for coalition in ACTIVE:
        assert (
            sum(projection[player] for player in range(N) if coalition >> player & 1)
            == game[coalition]
        )
    normals = tuple(normal(coalition) for coalition in ACTIVE)
    assert all(multiplier >= 0 for multiplier in multipliers)
    assert all(
        projection[player] - reference[player]
        == sum(
            multipliers[index] * normals[index][player]
            for index in range(len(ACTIVE))
        )
        for player in range(N)
    )


def main() -> int:
    verify_game(LOWER)
    verify_game(UPPER)
    assert all(
        UPPER[coalition] - LOWER[coalition]
        == (F(1, 1000) if coalition == CHANGED else 0)
        for coalition in range(GRAND + 1)
    )
    verify_projection(LOWER, LOWER_PROJECTION, MULTIPLIERS[0])
    verify_projection(UPPER, UPPER_PROJECTION, MULTIPLIERS[1])
    assert UPPER_PROJECTION[0] - LOWER_PROJECTION[0] == -F(1, 42000)
    print("PASS: exact games; exact KKT projections; protected drop = -1/42000")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
