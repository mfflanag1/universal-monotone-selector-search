#!/usr/bin/env python3
"""Verify and find the maximal domain of an affine splice-scale chart."""

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
Pair = tuple[F, F]


def members(coalition: int) -> tuple[int, ...]:
    return tuple(player for player in range(N) if coalition >> player & 1)


def add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def scale(value: Pair, coefficient: F) -> Pair:
    return value[0] * coefficient, value[1] * coefficient


def evaluate(value: Pair, parameter: F) -> F:
    return value[0] + value[1] * parameter


def affine(left: F, right: F, left_parameter: F, right_parameter: F) -> Pair:
    slope = (right - left) / (right_parameter - left_parameter)
    return left - slope * left_parameter, slope


def load(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    payload["_games"] = [
        tuple(F(value) for value in game) for game in payload["family"]["games"]
    ]
    payload["_edges"] = [
        (int(edge["lower"]), int(edge["upper"]), int(edge["coalition"]))
        for edge in payload["family"]["edges"]
    ]
    payload["_allocations"] = [
        tuple(F(value) for value in allocation)
        for allocation in payload["exact_allocations"]
    ]
    payload["_box_allocations"] = [
        tuple(F(value) for value in allocation)
        for allocation in payload["exact_box_allocations"]
    ]
    return payload


def common_budget_point(games: Sequence[Sequence[F]]) -> tuple[F, ...]:
    envelope = [max(game[coalition] for game in games) for coalition in range(32)]
    inequalities = [
        [F(-1) if coalition >> player & 1 else F(0) for player in range(N)]
        for coalition in range(1, GRAND)
    ]
    status, _, point = exact_linprog(
        [F(1)] * N,
        inequalities,
        [-envelope[coalition] for coalition in range(1, GRAND)],
    )
    if status != "optimal" or point is None:
        raise RuntimeError("common-budget LP failed")
    return tuple(point)


def dual_objective(
    certificate: dict[str, Any], games: Sequence[Sequence[Pair]]
) -> Pair:
    result = (F(0), F(0))
    for active in certificate["active_inequalities"]:
        row = active["row"]
        weight = F(active["weight"])
        if row[0] == "core":
            bound = scale(games[int(row[1])][int(row[2])], F(-1))
        elif row[0] == "box_lower":
            game, player = int(row[1]), int(row[2])
            bound = scale(games[game][1 << player], F(-1))
        elif row[0] == "box_upper":
            game, player = int(row[1]), int(row[2])
            bound = add(
                games[game][GRAND],
                scale(games[game][GRAND ^ (1 << player)], F(-1)),
            )
        elif row[0] == "monotonicity":
            bound = (F(0), F(0))
        else:
            raise RuntimeError(f"unknown dual row {row}")
        result = add(result, scale(bound, weight))
    for active in certificate["active_equalities"]:
        result = add(
            result,
            scale(games[int(active["game"])][GRAND], F(active["weight"])),
        )
    return result


def formula(value: Pair) -> str:
    constant, slope = value
    if slope == 0:
        return str(constant)
    return f"{constant}+({slope})*s"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("left_archive", type=Path)
    parser.add_argument("right_archive", type=Path)
    parser.add_argument("--left-parameter", type=F, required=True)
    parser.add_argument("--right-parameter", type=F, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not F(0) < args.left_parameter < args.right_parameter:
        raise ValueError("parameters must satisfy 0 < left < right")

    left = load(args.left_archive)
    right = load(args.right_archive)
    if left["_edges"] != right["_edges"]:
        raise RuntimeError("archive edge sets are not label-identical")
    if len(left["_games"]) != len(right["_games"]):
        raise RuntimeError("archive node counts differ")

    games = [
        tuple(
            affine(
                left["_games"][node][coalition],
                right["_games"][node][coalition],
                args.left_parameter,
                args.right_parameter,
            )
            for coalition in range(32)
        )
        for node in range(len(left["_games"]))
    ]
    allocations = [
        tuple(
            affine(
                left["_allocations"][node][player],
                right["_allocations"][node][player],
                args.left_parameter,
                args.right_parameter,
            )
            for player in range(N)
        )
        for node in range(len(games))
    ]
    box_allocations = [
        tuple(
            affine(
                left["_box_allocations"][node][player],
                right["_box_allocations"][node][player],
                args.left_parameter,
                args.right_parameter,
            )
            for player in range(N)
        )
        for node in range(len(games))
    ]
    margin = affine(
        F(left["max_min_monotonicity_margin_exact"]),
        F(right["max_min_monotonicity_margin_exact"]),
        args.left_parameter,
        args.right_parameter,
    )
    box_margin = affine(
        F(left["box_max_min_monotonicity_margin_exact"]),
        F(right["box_max_min_monotonicity_margin_exact"]),
        args.left_parameter,
        args.right_parameter,
    )
    left_common = common_budget_point(left["_games"])
    right_common = common_budget_point(right["_games"])
    common_point = tuple(
        affine(
            left_common[player],
            right_common[player],
            args.left_parameter,
            args.right_parameter,
        )
        for player in range(N)
    )
    common_budget = (F(0), F(0))
    for value in common_point:
        common_budget = add(common_budget, value)
    gap = add(common_budget, scale(games[0][GRAND], F(-1)))
    left_gap = F(left["common_core_budget_gap_exact"])
    right_gap = F(right["common_core_budget_gap_exact"])
    if gap != affine(
        left_gap,
        right_gap,
        args.left_parameter,
        args.right_parameter,
    ):
        raise RuntimeError("common-budget point does not interpolate the archived gaps")
    left_envelope = [
        max(game[coalition] for game in left["_games"])
        for coalition in range(32)
    ]
    complementary_pair: tuple[int, int, int, int] | None = None
    for coalition in range(1, GRAND):
        complement = GRAND ^ coalition
        if coalition >= complement:
            continue
        if (
            left_envelope[coalition]
            + left_envelope[complement]
            - left["_games"][0][GRAND]
            != left_gap
        ):
            continue
        first_node = next(
            node
            for node, game in enumerate(left["_games"])
            if game[coalition] == left_envelope[coalition]
        )
        second_node = next(
            node
            for node, game in enumerate(left["_games"])
            if game[complement] == left_envelope[complement]
        )
        lower_gap = add(
            add(games[first_node][coalition], games[second_node][complement]),
            scale(games[0][GRAND], F(-1)),
        )
        if lower_gap == gap:
            complementary_pair = (
                coalition,
                complement,
                first_node,
                second_node,
            )
            break
    if complementary_pair is None:
        raise RuntimeError("no fixed complementary lower certificate matches the gap")

    constraints: list[tuple[bool, str, Pair]] = []
    facets = load_facets()[0]
    for node, game in enumerate(games):
        if game[0] != (F(0), F(0)):
            raise RuntimeError("zero worth is not an affine identity")
        for coalition in range(32):
            for player in range(N):
                if coalition >> player & 1:
                    continue
                constraints.append(
                    (
                        False,
                        f"game_monotonicity:{node}:{coalition}:{player}",
                        add(
                            game[coalition | (1 << player)],
                            scale(game[coalition], F(-1)),
                        ),
                    )
                )
        for facet_index, facet in enumerate(facets):
            value = (F(0), F(0))
            for coalition, coefficient in enumerate(facet):
                if coefficient:
                    value = add(value, scale(game[coalition], F(coefficient)))
            constraints.append((False, f"exact_facet:{node}:{facet_index}", value))

    for edge_index, (lower, upper, coalition) in enumerate(left["_edges"]):
        delta = add(games[upper][coalition], scale(games[lower][coalition], F(-1)))
        constraints.append((True, f"edge:{edge_index}", delta))
        for target in range(32):
            if target == coalition:
                continue
            if add(games[upper][target], scale(games[lower][target], F(-1))) != (
                F(0),
                F(0),
            ):
                raise RuntimeError("edge changes a second worth parametrically")

    for node, (game, allocation, box_allocation) in enumerate(
        zip(games, allocations, box_allocations, strict=True)
    ):
        if sum((value[0] for value in allocation), F(0)) != game[GRAND][0] or sum(
            (value[1] for value in allocation), F(0)
        ) != game[GRAND][1]:
            raise RuntimeError("core efficiency is not an identity")
        if sum((value[0] for value in box_allocation), F(0)) != game[GRAND][0] or sum(
            (value[1] for value in box_allocation), F(0)
        ) != game[GRAND][1]:
            raise RuntimeError("box efficiency is not an identity")
        for coalition in range(1, GRAND):
            slack = scale(game[coalition], F(-1))
            for player in members(coalition):
                slack = add(slack, allocation[player])
            constraints.append((False, f"core:{node}:{coalition}", slack))
        for player in range(N):
            constraints.append(
                (
                    False,
                    f"box_lower:{node}:{player}",
                    add(box_allocation[player], scale(game[1 << player], F(-1))),
                )
            )
            constraints.append(
                (
                    False,
                    f"box_upper:{node}:{player}",
                    add(
                        add(game[GRAND], scale(game[GRAND ^ (1 << player)], F(-1))),
                        scale(box_allocation[player], F(-1)),
                    ),
                )
            )
        for coalition in range(1, GRAND):
            slack = scale(game[coalition], F(-1))
            for player in members(coalition):
                slack = add(slack, common_point[player])
            constraints.append((False, f"common_point:{node}:{coalition}", slack))

    for edge_index, (lower, upper, coalition) in enumerate(left["_edges"]):
        for player in members(coalition):
            core_gain = add(
                allocations[upper][player], scale(allocations[lower][player], F(-1))
            )
            constraints.append(
                (
                    False,
                    f"core_margin:{edge_index}:{player}",
                    add(core_gain, scale(margin, F(-1))),
                )
            )
            box_gain = add(
                box_allocations[upper][player],
                scale(box_allocations[lower][player], F(-1)),
            )
            constraints.append(
                (
                    False,
                    f"box_margin:{edge_index}:{player}",
                    add(box_gain, scale(box_margin, F(-1))),
                )
            )

    if dual_objective(left["exact_margin_dual_certificate"], games) != scale(
        margin, F(-1)
    ):
        raise RuntimeError("core dual objective is not affine-optimal")
    if dual_objective(left["exact_box_margin_dual_certificate"], games) != scale(
        box_margin, F(-1)
    ):
        raise RuntimeError("box dual objective is not affine-optimal")

    lower = F(0)
    upper: F | None = None
    lower_strict = True
    upper_strict = False
    lower_blockers = ["positive_scale"]
    upper_blockers: list[str] = []
    for strict, label, (constant, slope) in constraints:
        if slope == 0:
            if constant < 0 or (strict and constant == 0):
                raise RuntimeError(f"globally infeasible constraint {label}")
            continue
        root = -constant / slope
        if root < 0:
            continue
        if slope > 0:
            if root > lower:
                lower, lower_strict, lower_blockers = root, strict, [label]
            elif root == lower:
                lower_strict = lower_strict or strict
                lower_blockers.append(label)
        else:
            if upper is None or root < upper:
                upper, upper_strict, upper_blockers = root, strict, [label]
            elif root == upper:
                upper_strict = upper_strict or strict
                upper_blockers.append(label)
    if upper is not None and (
        lower > upper or (lower == upper and (lower_strict or upper_strict))
    ):
        raise RuntimeError("affine chart has no positive domain")
    for _, label, value in constraints:
        if evaluate(value, args.left_parameter) < 0 or evaluate(
            value, args.right_parameter
        ) < 0:
            raise RuntimeError(f"endpoint verification failed: {label}")

    result = {
        "status": "exact_affine_splice_scale_chart",
        "endpoint_parameters": [str(args.left_parameter), str(args.right_parameter)],
        "node_count": len(games),
        "edge_count": len(left["_edges"]),
        "constraint_count": len(constraints),
        "formulas": {
            "core_margin": formula(margin),
            "box_margin": formula(box_margin),
            "common_core_gap": formula(gap),
        },
        "common_core_lower_certificate": {
            "coalitions": [complementary_pair[0], complementary_pair[1]],
            "nodes": [complementary_pair[2], complementary_pair[3]],
            "weights": ["1", "1"],
        },
        "maximal_certificate_domain": {
            "lower": str(lower),
            "lower_strict": lower_strict,
            "upper": None if upper is None else str(upper),
            "upper_strict": upper_strict,
            "lower_blockers": lower_blockers,
            "upper_blockers": upper_blockers,
        },
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
