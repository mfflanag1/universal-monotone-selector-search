#!/usr/bin/env python3
"""Optimize a transverse/template splice against a prescribed dual mixture."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from family_search import box_family_slack_float, common_core_gap_float
from n5_facet_search import (
    ExactConeModel,
    flow_objective,
    margin_dual,
    permute_mask,
    solve_model,
)


F = Fraction
GRAND = 31


def record_from_archive(payload: dict) -> dict:
    if payload.get("record") is not None:
        return payload["record"]
    if payload.get("search") is not None:
        return payload["search"]
    if payload.get("family") is not None:
        family = payload["family"]
        normalizer = F(family["games"][0][GRAND])
        return {
            "edges": [
                [
                    int(edge["lower"]),
                    int(edge["upper"]),
                    int(edge["coalition"]),
                ]
                for edge in family["edges"]
            ],
            "bumps": [
                float(F(edge["delta"]) / normalizer)
                for edge in family["edges"]
            ],
        }
    return payload


def add_dual_objective(
    objective: np.ndarray,
    model: ExactConeModel,
    archive: dict,
    coefficient: float,
    node_map: dict[int, int] | None = None,
    permutation: tuple[int, ...] | None = None,
) -> None:
    dual = archive["exact_margin_dual_certificate"]
    for active in dual["active_inequalities"]:
        row = active["row"]
        if row[0] != "core":
            continue
        node = int(row[1])
        coalition = int(row[2])
        if node_map is not None:
            node = node_map[node]
        if permutation is not None:
            coalition = permute_mask(coalition, permutation)
        objective[model.game_index(node, coalition)] += (
            coefficient * float(F(active["weight"]))
        )
    for active in dual["active_equalities"]:
        node = int(active["game"])
        if node_map is not None:
            node = node_map[node]
        objective[model.game_index(node, GRAND)] -= (
            coefficient * float(F(active["weight"]))
        )


def fixed_dual_value(
    games: list[list[float]],
    archive: dict,
    coefficient: float,
    node_map: dict[int, int] | None = None,
    permutation: tuple[int, ...] | None = None,
) -> float:
    value = 0.0
    dual = archive["exact_margin_dual_certificate"]
    for active in dual["active_inequalities"]:
        row = active["row"]
        if row[0] != "core":
            continue
        node = int(row[1])
        coalition = int(row[2])
        if node_map is not None:
            node = node_map[node]
        if permutation is not None:
            coalition = permute_mask(coalition, permutation)
        value += (
            coefficient
            * float(F(active["weight"]))
            * games[node][coalition]
        )
    for active in dual["active_equalities"]:
        node = int(active["game"])
        if node_map is not None:
            node = node_map[node]
        value -= (
            coefficient
            * float(F(active["weight"]))
            * games[node][GRAND]
        )
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attachment", type=int, required=True)
    parser.add_argument("--template-anchor", type=int, required=True)
    parser.add_argument("--permutation", default="2,1,4,0,3")
    parser.add_argument("--scale", type=float, default=0.00005)
    parser.add_argument(
        "--alphas", default="0.01,0.03,0.1,0.3,1,3,10,30,100"
    )
    parser.add_argument("--iterations", type=int, default=0)
    parser.add_argument("--envelope-nodes")
    parser.add_argument("--envelope-coalitions")
    parser.add_argument("--envelope-gap", type=float, default=0.0)
    parser.add_argument("--envelope-multiplier", type=float, default=1.0)
    args = parser.parse_args()

    permutation = tuple(
        int(value) for value in args.permutation.split(",")
    )
    if sorted(permutation) != list(range(5)):
        raise ValueError("permutation must contain 0,1,2,3,4 exactly once")
    alphas = [float(value) for value in args.alphas.split(",")]

    parent_archive = json.loads(args.parent.read_text())
    parent = record_from_archive(parent_archive)
    template_archive = json.loads(args.template.read_text())
    template = template_archive["family"]
    edges = [tuple(int(value) for value in edge) for edge in parent["edges"]]
    bumps = [float(value) for value in parent["bumps"]]
    game_count = 1 + max(
        max(lower, upper) for lower, upper, _ in edges
    )
    if not 0 <= args.attachment < game_count:
        raise ValueError("attachment is outside the parent game range")

    node_map = {args.template_anchor: args.attachment}
    next_node = game_count
    for node in range(len(template["games"])):
        if node == args.template_anchor:
            continue
        node_map[node] = next_node
        next_node += 1
    for edge in template["edges"]:
        edges.append(
            (
                node_map[int(edge["lower"])],
                node_map[int(edge["upper"])],
                permute_mask(int(edge["coalition"]), permutation),
            )
        )
        bumps.append(args.scale * float(F(edge["delta"])))

    model = ExactConeModel(edges, bumps)
    for game in range(1, model.game_count):
        model.add_eq({model.game_index(game, GRAND): 1.0}, 1.0)
    model.a_eq = model.sparse(model.eq_rows)
    if args.envelope_nodes is not None:
        envelope_nodes = [
            int(value) for value in args.envelope_nodes.split(",")
        ]
        if args.envelope_coalitions is None:
            envelope_coalitions = [1 << player for player in range(5)]
        else:
            envelope_coalitions = [
                int(value)
                for value in args.envelope_coalitions.split(",")
            ]
        if len(envelope_nodes) != len(envelope_coalitions):
            raise ValueError("envelope node and coalition counts differ")
        row: dict[int, float] = {}
        for node, coalition in zip(
            envelope_nodes, envelope_coalitions, strict=True
        ):
            column = model.game_index(node, coalition)
            row[column] = row.get(column, 0.0) - 1.0
        model.add_ub(
            row,
            -args.envelope_multiplier * (1.0 + args.envelope_gap),
        )
        model.a_ub = model.sparse(model.ub_rows)
    rows = []
    candidates = []
    for alpha in alphas:
        objective = np.zeros(model.variable_count)
        add_dual_objective(objective, model, parent_archive, 1.0)
        add_dual_objective(
            objective,
            model,
            template_archive,
            alpha,
            node_map,
            permutation,
        )
        games = solve_model(model, objective)
        fixed_bound = (
            fixed_dual_value(games, parent_archive, 1.0)
            + fixed_dual_value(
                games,
                template_archive,
                alpha,
                node_map,
                permutation,
            )
        ) / (1.0 + alpha)
        row = {"alpha": alpha, "fixed_dual_margin_bound": fixed_bound}
        rows.append(row)
        candidates.append((fixed_bound, alpha, games))
        print(json.dumps(row), flush=True)

    fixed_bound, alpha, games = min(candidates)
    common_grand = games[0][GRAND]
    grand_residual = max(
        abs(game[GRAND] - common_grand) for game in games
    )
    if grand_residual > 1e-6:
        raise RuntimeError(
            "glued topology does not preserve common grand: "
            f"max residual {grand_residual}"
        )
    for game in games:
        game[GRAND] = common_grand
    margin, iw, ew, metadata = margin_dual(games, edges)
    alternating_trace = [margin]
    for _ in range(args.iterations):
        objective = flow_objective(model, iw, ew, metadata)
        next_games = solve_model(model, objective)
        for game in next_games:
            game[GRAND] = common_grand
        next_margin, next_iw, next_ew, next_metadata = margin_dual(
            next_games, edges
        )
        alternating_trace.append(next_margin)
        if next_margin >= margin - 1e-10:
            break
        games = next_games
        margin = next_margin
        iw, ew, metadata = next_iw, next_ew, next_metadata
    box_margin, _ = box_family_slack_float(games, edges, 5)
    result = {
        "status": "targeted_dual_splice_float",
        "parent": str(args.parent),
        "template": str(args.template),
        "attachment": args.attachment,
        "template_anchor": args.template_anchor,
        "permutation": list(permutation),
        "scale": args.scale,
        "rows": rows,
        "best_alpha": alpha,
        "best_fixed_dual_margin_bound": fixed_bound,
        "best_margin_float": margin,
        "alternating_margin_trace": alternating_trace,
        "box_margin_float": box_margin,
        "non_atomic_facet_tax_float": (
            box_margin - margin if box_margin is not None else None
        ),
        "common_core_gap_float": common_core_gap_float(games, 5),
        "envelope": {
            "nodes": args.envelope_nodes,
            "coalitions": args.envelope_coalitions,
            "gap": args.envelope_gap,
            "multiplier": args.envelope_multiplier,
        },
        "record": {
            "n": 5,
            "edges": [list(edge) for edge in edges],
            "bumps": bumps,
            "best_margin_float": margin,
            "best_games_float": games,
            "facet_count": len(model.facets),
            "variables": model.variable_count,
            "inequalities": len(model.ub_rows),
            "equalities": len(model.eq_rows),
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "best_alpha",
                    "best_fixed_dual_margin_bound",
                    "best_margin_float",
                    "box_margin_float",
                    "non_atomic_facet_tax_float",
                    "common_core_gap_float",
                )
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
