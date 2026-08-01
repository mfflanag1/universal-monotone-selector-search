#!/usr/bin/env python3
"""Find the maximal positive-parameter domain of a rolling-triple chart."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import n5_parametric_rolling_triple_verify as chart
from n5_facet_search import load_facets


F = Fraction
Pair = tuple[F, F]
N = 5
GRAND = 31


def add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def scale(value: Pair, coefficient: F) -> Pair:
    return coefficient * value[0], coefficient * value[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("left_archive", type=Path)
    parser.add_argument("right_archive", type=Path)
    parser.add_argument("--left-parameter", type=F, required=True)
    parser.add_argument("--right-parameter", type=F, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 0 < args.left_parameter < args.right_parameter:
        raise ValueError("parameters must satisfy 0 < left < right")

    chart.LEFT_PARAMETER = args.left_parameter
    chart.RIGHT_PARAMETER = args.right_parameter
    first = chart.load(args.left_archive)
    second = chart.load(args.right_archive)
    mapping = chart.candidate_mapping(first, second)
    facets = load_facets()[0]

    games = [
        tuple(
            chart.affine(
                first["_games"][node][coalition],
                second["_games"][mapping[node]][coalition],
            )
            for coalition in range(1 << N)
        )
        for node in range(len(first["_games"]))
    ]
    allocations = [
        tuple(
            chart.affine(
                first["_allocations"][node][player],
                second["_allocations"][mapping[node]][player],
            )
            for player in range(N)
        )
        for node in range(len(first["_games"]))
    ]

    constraints: list[tuple[bool, str, Pair]] = []
    for node, game in enumerate(games):
        for coalition in range(1 << N):
            for player in range(N):
                if coalition >> player & 1:
                    continue
                upper = coalition | (1 << player)
                constraints.append(
                    (
                        False,
                        f"monotonicity:{node}:{coalition}:{upper}",
                        add(game[upper], scale(game[coalition], F(-1))),
                    )
                )
        for facet_index, facet in enumerate(facets):
            value = (F(0), F(0))
            for coalition in range(1, 1 << N):
                if facet[coalition]:
                    value = add(
                        value,
                        scale(game[coalition], F(facet[coalition])),
                    )
            constraints.append(
                (False, f"exact_facet:{node}:{facet_index}", value)
            )

    for edge_index, (lower, upper, coalition) in enumerate(first["_edges"]):
        constraints.append(
            (
                True,
                f"edge:{edge_index}:{lower}:{upper}:{coalition}",
                add(
                    games[upper][coalition],
                    scale(games[lower][coalition], F(-1)),
                ),
            )
        )

    margin = (F(0), F(2, 33))
    for node, (game, point) in enumerate(
        zip(games, allocations, strict=True)
    ):
        efficiency = (F(0), F(0))
        for value in point:
            efficiency = add(efficiency, value)
        if efficiency != game[GRAND]:
            raise RuntimeError(f"parametric efficiency identity failed: {node}")
        for coalition in range(1, GRAND):
            slack = scale(game[coalition], F(-1))
            for player in range(N):
                if coalition >> player & 1:
                    slack = add(slack, point[player])
            constraints.append(
                (False, f"core:{node}:{coalition}", slack)
            )

    for edge_index, (lower, upper, coalition) in enumerate(first["_edges"]):
        for player in range(N):
            if coalition >> player & 1:
                gain = add(
                    allocations[upper][player],
                    scale(allocations[lower][player], F(-1)),
                )
                constraints.append(
                    (
                        False,
                        f"margin:{edge_index}:{player}",
                        add(gain, scale(margin, F(-1))),
                    )
                )

    common_point = tuple(map(F, (1, 1, 3, 3, 0)))
    for node, game in enumerate(games):
        for coalition in range(1, GRAND):
            value = sum(
                (
                    common_point[player]
                    for player in range(N)
                    if coalition >> player & 1
                ),
                F(0),
            )
            constraints.append(
                (
                    False,
                    f"common_point:{node}:{coalition}",
                    add(scale(game[coalition], F(-1)), (value, F(0))),
                )
            )

    lower_z = F(0)
    upper_z: F | None = None
    lower_strict = True
    upper_strict = False
    lower_labels: list[str] = ["positive_parameter"]
    upper_labels: list[str] = []
    for strict, label, (constant, inverse) in constraints:
        if inverse == 0:
            if constant < 0 or (strict and constant == 0):
                raise RuntimeError(f"globally infeasible constraint: {label}")
            continue
        root = -constant / inverse
        if inverse > 0 and root >= 0:
            if root > lower_z:
                lower_z = root
                lower_strict = strict
                lower_labels = [label]
            elif root == lower_z:
                lower_strict = lower_strict or strict
                lower_labels.append(label)
        elif inverse < 0 and root >= 0:
            if upper_z is None or root < upper_z:
                upper_z = root
                upper_strict = strict
                upper_labels = [label]
            elif root == upper_z:
                upper_strict = upper_strict or strict
                upper_labels.append(label)

    if upper_z is not None and (
        lower_z > upper_z
        or (lower_z == upper_z and (lower_strict or upper_strict))
    ):
        raise RuntimeError("the chart has no positive-parameter domain")

    lower_l = F(0) if upper_z is None else 1 / upper_z
    upper_l: F | None = None if lower_z == 0 else 1 / lower_z
    result = {
        "status": "maximal_parametric_chart_domain",
        "endpoint_parameters": [
            str(args.left_parameter),
            str(args.right_parameter),
        ],
        "inverse_parameter_domain": {
            "lower": str(lower_z),
            "lower_strict": lower_strict,
            "upper": None if upper_z is None else str(upper_z),
            "upper_strict": upper_strict,
        },
        "parameter_domain": {
            "lower": str(lower_l),
            "lower_strict": upper_strict,
            "upper": None if upper_l is None else str(upper_l),
            "upper_strict": lower_strict,
        },
        "lower_parameter_blockers": upper_labels,
        "upper_parameter_blockers": lower_labels,
        "constraint_count": len(constraints),
        "node_count": len(games),
        "edge_count": len(first["_edges"]),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
