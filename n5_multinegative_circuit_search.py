#!/usr/bin/env python3
"""Optimize protected routings induced by a non-semi-balanced incidence circuit."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from n5_mixed_terminal_search import GameCone, Topology


N = 5
GRAND = (1 << N) - 1


def player_occurrences(blocks: tuple[int, ...], player: int) -> tuple[int, ...]:
    return tuple(index for index, block in enumerate(blocks) if block >> player & 1)


def routings(
    positive: tuple[int, ...], negative: tuple[int, ...]
):
    choices = []
    for player in range(N):
        sources = player_occurrences(positive, player)
        sinks = player_occurrences(negative, player)
        if len(sources) != len(sinks):
            raise ValueError(f"incidence is not balanced for player {player}")
        choices.append(
            tuple(
                tuple(zip(sources, permutation, strict=True))
                for permutation in itertools.permutations(sinks)
            )
        )
    seen = set()
    for player_matchings in itertools.product(*choices):
        edge_players: dict[tuple[int, int], int] = {}
        for player, matching in enumerate(player_matchings):
            for source, sink in matching:
                pair = (source, len(positive) + sink)
                edge_players[pair] = edge_players.get(pair, 0) | (1 << player)
        arcs = tuple(
            (source, sink, players)
            for (source, sink), players in sorted(edge_players.items())
        )
        if arcs in seen:
            continue
        seen.add(arcs)
        yield arcs


def topology(
    positive: tuple[int, ...], negative: tuple[int, ...], arcs: tuple[tuple[int, int, int], ...]
) -> Topology:
    divergence = [[0] * N for _ in range(len(positive) + len(negative))]
    for source, sink, players in arcs:
        for player in range(N):
            if players >> player & 1:
                divergence[source][player] += 1
                divergence[sink][player] -= 1
    expected = tuple(
        tuple(int(block >> player & 1) for player in range(N))
        for block in positive
    ) + tuple(
        tuple(-int(block >> player & 1) for player in range(N))
        for block in negative
    )
    if tuple(tuple(row) for row in divergence) != expected:
        raise RuntimeError("routing does not realize the requested terminal divergences")
    return Topology(len(expected), arcs, expected)


def objective(cone: GameCone, positive: tuple[int, ...], negative: tuple[int, ...]) -> np.ndarray:
    coefficients = np.zeros(cone.variable_count)
    for node, block in enumerate(positive):
        coefficients[cone.column(node, block)] += 1.0
    offset = len(positive)
    for index, block in enumerate(negative):
        complement = GRAND ^ block
        if complement:
            coefficients[cone.column(offset + index, complement)] += 1.0
        coefficients[cone.column(offset + index, GRAND)] -= 1.0
    return coefficients


def solve(one_routing: Topology, positive: tuple[int, ...], negative: tuple[int, ...]):
    cone = GameCone(N, one_routing)
    coefficients = objective(cone, positive, negative)
    result = linprog(
        -coefficients,
        A_ub=cone.a_ub,
        b_ub=cone.b_ub,
        A_eq=cone.a_eq,
        b_eq=cone.b_eq,
        bounds=[(0.0, 1.0)] * cone.variable_count,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return -float(result.fun), result.x


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    # Example 3 of Studeny--Kratochvil has the irreducible two-negative
    # dependence -ab-ac+ad+abc+abce=N.
    positive = (0b01001, 0b00111, 0b10111)  # ad, abc, abce
    negative = (0b00011, 0b00101, GRAND)  # ab, ac, N
    records = []
    for index, arcs in enumerate(routings(positive, negative)):
        one_routing = topology(positive, negative, arcs)
        optimum, point = solve(one_routing, positive, negative)
        record = {
            "routing": index,
            "arcs": [list(arc) for arc in arcs],
            "optimum": optimum,
        }
        if optimum > 1e-9:
            record["games"] = [
                [
                    float(value)
                    for value in point[node * (GRAND + 1) : (node + 1) * (GRAND + 1)]
                ]
                for node in range(one_routing.node_count)
            ]
        records.append(record)
        print(json.dumps(record | {"games": "omitted"} if "games" in record else record), flush=True)

    best = max(records, key=lambda row: row["optimum"])
    payload = {
        "status": "positive_found" if best["optimum"] > 1e-9 else "globally_nonpositive",
        "positive_blocks": list(positive),
        "negative_blocks": list(negative),
        "routing_count": len(records),
        "best": best,
        "records": [{key: value for key, value in row.items() if key != "games"} for row in records],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "routing_count", "best") if key in payload}, indent=2))
    return 1 if payload["status"] == "positive_found" else 0


if __name__ == "__main__":
    raise SystemExit(main())
