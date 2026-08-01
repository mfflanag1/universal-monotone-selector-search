#!/usr/bin/env python3
"""Independently verify an archived sparse exact n=5 family certificate."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

from exact_lp import exact_linprog
from n5_facet_search import load_facets


F = Fraction
N = 5
GRAND = (1 << N) - 1


def members(coalition: int) -> tuple[int, ...]:
    return tuple(player for player in range(N) if coalition >> player & 1)


def common_core_gap(games: Sequence[Sequence[F]]) -> F:
    envelope = [max(game[coalition] for game in games) for coalition in range(32)]
    inequalities = [
        [F(-1) if coalition >> player & 1 else F(0) for player in range(N)]
        for coalition in range(1, GRAND)
    ]
    status, optimum, _ = exact_linprog(
        [F(1)] * N,
        inequalities,
        [-envelope[coalition] for coalition in range(1, GRAND)],
    )
    if status != "optimal" or optimum is None:
        raise AssertionError("common-core LP did not solve")
    return optimum - games[0][GRAND]


def row_entries(metadata: Sequence[Any], game_count: int) -> tuple[dict[int, F], F]:
    tag = metadata[0]
    margin_index = game_count * N
    if tag == "core":
        game, coalition = int(metadata[1]), int(metadata[2])
        return (
            {game * N + player: F(-1) for player in members(coalition)},
            None,
        )
    if tag == "box_lower":
        game, player = int(metadata[1]), int(metadata[2])
        return ({game * N + player: F(-1)}, None)
    if tag == "box_upper":
        game, player = int(metadata[1]), int(metadata[2])
        return ({game * N + player: F(1)}, None)
    if tag == "monotonicity":
        lower, upper, player = (
            int(metadata[1]),
            int(metadata[2]),
            int(metadata[4]),
        )
        return (
            {
                lower * N + player: F(1),
                upper * N + player: F(-1),
                margin_index: F(1),
            },
            F(0),
        )
    raise AssertionError(f"unknown inequality row {metadata}")


def row_rhs(
    metadata: Sequence[Any], games: Sequence[Sequence[F]]
) -> F:
    tag = metadata[0]
    if tag == "core":
        return -games[int(metadata[1])][int(metadata[2])]
    if tag == "box_lower":
        game, player = int(metadata[1]), int(metadata[2])
        return -games[game][1 << player]
    if tag == "box_upper":
        game, player = int(metadata[1]), int(metadata[2])
        return games[game][GRAND] - games[game][GRAND ^ (1 << player)]
    if tag == "monotonicity":
        return F(0)
    raise AssertionError(f"unknown inequality row {metadata}")


def verify_allocations(
    games: Sequence[Sequence[F]],
    edges: Sequence[tuple[int, int, int]],
    allocations: Sequence[Sequence[F]],
    margin: F,
    box: bool,
) -> None:
    assert len(allocations) == len(games)
    for game, allocation in zip(games, allocations, strict=True):
        assert sum(allocation, F(0)) == game[GRAND]
        if box:
            for player in range(N):
                assert allocation[player] >= game[1 << player]
                assert allocation[player] <= (
                    game[GRAND] - game[GRAND ^ (1 << player)]
                )
        else:
            for coalition in range(1, GRAND):
                assert sum(
                    (allocation[player] for player in members(coalition)), F(0)
                ) >= game[coalition]
    for lower, upper, coalition in edges:
        for player in members(coalition):
            assert allocations[upper][player] >= allocations[lower][player] + margin


def verify_dual(
    games: Sequence[Sequence[F]],
    certificate: dict[str, Any],
    margin: F,
) -> None:
    variable_count = len(games) * N + 1
    stationarity = [F(0)] * variable_count
    dual_objective = F(0)
    for active in certificate["active_inequalities"]:
        weight = F(active["weight"])
        assert weight <= 0
        entries, _ = row_entries(active["row"], len(games))
        for column, coefficient in entries.items():
            stationarity[column] += weight * coefficient
        dual_objective += weight * row_rhs(active["row"], games)
    for active in certificate["active_equalities"]:
        game = int(active["game"])
        weight = F(active["weight"])
        for player in range(N):
            stationarity[game * N + player] += weight
        dual_objective += weight * games[game][GRAND]
    objective = [F(0)] * variable_count
    objective[-1] = F(-1)
    assert stationarity == objective
    assert dual_objective == -margin
    assert F(certificate["objective"]) == -margin


def verify_connected(game_count: int, edges: Sequence[tuple[int, int, int]]) -> None:
    adjacency = [set() for _ in range(game_count)]
    for lower, upper, _ in edges:
        adjacency[lower].add(upper)
        adjacency[upper].add(lower)
    reached = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for neighbor in adjacency[node] - reached:
            reached.add(neighbor)
            stack.append(neighbor)
    assert len(reached) == game_count


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.archive.read_text())
    assert int(payload["n"]) == N
    games = [tuple(F(value) for value in game) for game in payload["family"]["games"]]
    edges = [
        (int(edge["lower"]), int(edge["upper"]), int(edge["coalition"]))
        for edge in payload["family"]["edges"]
    ]
    assert len(games) == int(payload["node_count"])
    assert len(edges) == int(payload["edge_count"])
    assert len({game[GRAND] for game in games}) == 1
    for lower, upper, coalition in edges:
        differences = [
            games[upper][target] - games[lower][target]
            for target in range(GRAND + 1)
        ]
        assert differences[coalition] > 0
        assert all(
            value == 0
            for target, value in enumerate(differences)
            if target != coalition
        )
    facets, _ = load_facets()
    for game in games:
        assert all(
            game[coalition] <= game[coalition | (1 << player)]
            for coalition in range(GRAND + 1)
            for player in range(N)
            if not coalition >> player & 1
        )
        assert all(
            sum(
                (F(coefficient) * game[index] for index, coefficient in enumerate(facet)),
                F(0),
            )
            >= 0
            for facet in facets
        )
    verify_connected(len(games), edges)
    gap = common_core_gap(games)
    assert gap == F(payload["common_core_budget_gap_exact"])
    margin = F(payload["max_min_monotonicity_margin_exact"])
    allocations = [
        tuple(F(value) for value in allocation)
        for allocation in payload["exact_allocations"]
    ]
    verify_allocations(games, edges, allocations, margin, False)
    verify_dual(games, payload["exact_margin_dual_certificate"], margin)
    box_margin = F(payload["box_max_min_monotonicity_margin_exact"])
    box_allocations = [
        tuple(F(value) for value in allocation)
        for allocation in payload["exact_box_allocations"]
    ]
    verify_allocations(games, edges, box_allocations, box_margin, True)
    verify_dual(games, payload["exact_box_margin_dual_certificate"], box_margin)
    assert box_margin - margin == F(payload["non_atomic_facet_tax_exact"])
    print(
        "PASS: "
        f"{len(games)} games, {len(edges)} edges, one component, "
        f"gap={gap}, margin={margin}, box={box_margin}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
