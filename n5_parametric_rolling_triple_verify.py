#!/usr/bin/env python3
"""Verify an all-parameter family interpolated from two rolling triples."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import networkx as nx
from scipy.optimize import linprog

from n5_facet_search import load_facets


F = Fraction
Pair = tuple[F, F]
N = 5
GRAND = 31
LEFT_PARAMETER = F(5)
RIGHT_PARAMETER = F(6)


def add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def scale(value: Pair, coefficient: F) -> Pair:
    return coefficient * value[0], coefficient * value[1]


def affine(first: F, second: F) -> Pair:
    """Return a+b/L through the two configured endpoints."""

    b = (first - second) / (
        1 / LEFT_PARAMETER - 1 / RIGHT_PARAMETER
    )
    return first - b / LEFT_PARAMETER, b


def nonnegative(value: Pair) -> bool:
    """Check a+b/L >= 0 throughout the configured interval."""

    return (
        value[0] + value[1] / LEFT_PARAMETER >= 0
        and value[0] + value[1] / RIGHT_PARAMETER >= 0
    )


def positive(value: Pair) -> bool:
    """Check a+b/L > 0 throughout the configured interval."""

    return (
        value[0] + value[1] / LEFT_PARAMETER > 0
        and value[0] + value[1] / RIGHT_PARAMETER > 0
    )


def load(path: Path) -> dict:
    payload = json.loads(path.read_text())
    payload["_games"] = [
        tuple(F(value) for value in game)
        for game in payload["family"]["games"]
    ]
    payload["_edges"] = [
        (
            int(edge["lower"]),
            int(edge["upper"]),
            int(edge["coalition"]),
        )
        for edge in payload["family"]["edges"]
    ]
    payload["_allocations"] = [
        tuple(F(value) for value in point)
        for point in payload["exact_allocations"]
    ]
    return payload


def graph(payload: dict) -> nx.DiGraph:
    result = nx.DiGraph()
    result.add_nodes_from(
        (node, {"index": node})
        for node in range(len(payload["_games"]))
    )
    for lower, upper, coalition in payload["_edges"]:
        result.add_edge(lower, upper, coalition=coalition)
    return result


def candidate_mapping(
    first: dict,
    second: dict,
) -> dict[int, int]:
    matcher = nx.algorithms.isomorphism.DiGraphMatcher(
        graph(first),
        graph(second),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "coalition", None
        ),
    )
    for mapping in matcher.isomorphisms_iter():
        if mapping.get(0) == 0:
            return mapping
    raise RuntimeError("rolling-triple graphs are not label-isomorphic")


def main() -> int:
    global LEFT_PARAMETER, RIGHT_PARAMETER

    parser = argparse.ArgumentParser()
    parser.add_argument("left_archive", type=Path)
    parser.add_argument("right_archive", type=Path)
    parser.add_argument(
        "--left-parameter", type=F, default=F(5)
    )
    parser.add_argument(
        "--right-parameter", type=F, default=F(6)
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 0 < args.left_parameter < args.right_parameter:
        raise ValueError("parameters must satisfy 0 < left < right")
    LEFT_PARAMETER = args.left_parameter
    RIGHT_PARAMETER = args.right_parameter

    first = load(args.left_archive)
    second = load(args.right_archive)
    facets = load_facets()[0]
    mapping = candidate_mapping(first, second)
    games = [
        tuple(
            affine(
                first["_games"][node][coalition],
                second["_games"][mapping[node]][coalition],
            )
            for coalition in range(1 << N)
        )
        for node in range(len(first["_games"]))
    ]
    allocations = [
        tuple(
            affine(
                first["_allocations"][node][player],
                second["_allocations"][mapping[node]][player],
            )
            for player in range(N)
        )
        for node in range(len(first["_games"]))
    ]

    for game in games:
        if game[0] != (F(0), F(0)):
            raise RuntimeError("zero normalization is not parametric")
        for coalition in range(1 << N):
            for player in range(N):
                if coalition >> player & 1:
                    continue
                upper = coalition | (1 << player)
                difference = add(
                    game[upper], scale(game[coalition], F(-1))
                )
                if not nonnegative(difference):
                    raise RuntimeError("parametric monotonicity failed")
        for facet in facets:
            value = (F(0), F(0))
            for coalition in range(1, 1 << N):
                if facet[coalition]:
                    value = add(
                        value,
                        scale(
                            game[coalition],
                            F(facet[coalition]),
                        ),
                    )
            if not nonnegative(value):
                raise RuntimeError("parametric exact facet failed")

    for lower, upper, coalition in first["_edges"]:
        delta = add(
            games[upper][coalition],
            scale(games[lower][coalition], F(-1)),
        )
        if not positive(delta):
            raise RuntimeError("parametric edge is not strictly positive")
        for target in range(1 << N):
            if target == coalition:
                continue
            difference = add(
                games[upper][target],
                scale(games[lower][target], F(-1)),
            )
            if difference != (F(0), F(0)):
                raise RuntimeError("parametric edge changes two worths")

    margin = (F(0), F(2, 33))
    for node, (game, point) in enumerate(
        zip(games, allocations, strict=True)
    ):
        efficiency = (F(0), F(0))
        for value in point:
            efficiency = add(efficiency, value)
        if efficiency != game[GRAND]:
            raise RuntimeError(f"efficiency failed at node {node}")
        for coalition in range(1, GRAND):
            slack = scale(game[coalition], F(-1))
            for player in range(N):
                if coalition >> player & 1:
                    slack = add(slack, point[player])
            if not nonnegative(slack):
                raise RuntimeError("parametric core allocation failed")
    for lower, upper, coalition in first["_edges"]:
        for player in range(N):
            if not coalition >> player & 1:
                continue
            gain = add(
                allocations[upper][player],
                scale(allocations[lower][player], F(-1)),
            )
            slack = add(gain, scale(margin, F(-1)))
            if not nonnegative(slack):
                raise RuntimeError("parametric margin primal failed")

    dual = first["exact_margin_dual_certificate"]
    dual_objective = (F(0), F(0))
    for active in dual["active_inequalities"]:
        row = active["row"]
        weight = F(active["weight"])
        if row[0] == "core":
            bound = scale(
                games[int(row[1])][int(row[2])], F(-1)
            )
        elif row[0] == "monotonicity":
            bound = (F(0), F(0))
        else:
            raise RuntimeError("unexpected core-dual row")
        dual_objective = add(
            dual_objective, scale(bound, weight)
        )
    for active in dual["active_equalities"]:
        node = int(active["game"])
        weight = F(active["weight"])
        dual_objective = add(
            dual_objective, scale(games[node][GRAND], weight)
        )
    if dual_objective != scale(margin, F(-1)):
        raise RuntimeError("parametric dual objective failed")

    common_point = tuple(map(F, (1, 1, 3, 3, 0)))
    for game in games:
        for coalition in range(1, GRAND):
            slack = scale(game[coalition], F(-1))
            value = sum(
                (
                    common_point[player]
                    for player in range(N)
                    if coalition >> player & 1
                ),
                F(0),
            )
            slack = add(slack, (value, F(0)))
            if not nonnegative(slack):
                raise RuntimeError("common budget-8 point failed")

    first_games = first["_games"]
    envelope = [
        max(game[coalition] for game in first_games)
        for coalition in range(1 << N)
    ]
    envelope_nodes = [
        next(
            node
            for node, game in enumerate(first_games)
            if game[coalition] == envelope[coalition]
        )
        for coalition in range(1 << N)
    ]
    common_result = linprog(
        [1.0] * N,
        A_ub=[
            [
                -1.0 if coalition >> player & 1 else 0.0
                for player in range(N)
            ]
            for coalition in range(1, GRAND)
        ],
        b_ub=[
            -float(envelope[coalition])
            for coalition in range(1, GRAND)
        ],
        bounds=[(None, None)] * N,
        method="highs",
    )
    if not common_result.success:
        raise RuntimeError("common-budget dual solve failed")
    common_weights = [
        -F(value).limit_denominator(1_000_000)
        for value in common_result.ineqlin.marginals
    ]
    for player in range(N):
        coverage = sum(
            (
                common_weights[coalition - 1]
                for coalition in range(1, GRAND)
                if coalition >> player & 1
            ),
            F(0),
        )
        if coverage != 1:
            raise RuntimeError("common-budget dual is not balanced")
    common_lower = (F(0), F(0))
    for coalition in range(1, GRAND):
        weight = common_weights[coalition - 1]
        if not weight:
            continue
        common_lower = add(
            common_lower,
            scale(
                games[envelope_nodes[coalition]][coalition],
                weight,
            ),
        )
    if common_lower != (F(8), F(0)):
        raise RuntimeError("common-budget lower certificate failed")

    result = {
        "status": "parametric_rolling_triple_verified",
        "parameter_domain": (
            f"real {LEFT_PARAMETER} <= L <= {RIGHT_PARAMETER}"
        ),
        "node_count": len(games),
        "edge_count": len(first["_edges"]),
        "exact_margin": "2/(33 L)",
        "grand_worth": "6 + 4/L",
        "common_core_upper_budget": "8",
        "common_core_gap": "2 - 4/L",
        "mapping": [mapping[node] for node in range(len(games))],
        "checks": {
            "all_280_exact_facets": True,
            "game_monotonicity": True,
            "strict_single_coordinate_edges": True,
            "parametric_core_primal": True,
            "parametric_core_dual": True,
            "common_budget_8_point": True,
            "common_budget_8_lower_certificate": True,
        },
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
