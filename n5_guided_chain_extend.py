#!/usr/bin/env python3
"""Extend a float family and seed optimization with its active margin dual."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from family_search import common_core_gap_float
from n5_facet_search import (
    ExactConeModel,
    flow_objective,
    margin_dual,
    permute_mask,
    solve_model,
)


F = Fraction


def load_payload(path: Path) -> tuple[dict, dict]:
    payload = json.loads(path.read_text())
    record = payload.get("record", payload)
    if record.get("best_games_float") is None and payload.get("family"):
        family = payload["family"]
        record = {
            "edges": [
                [
                    int(edge["lower"]),
                    int(edge["upper"]),
                    int(edge["coalition"]),
                ]
                for edge in family["edges"]
            ],
            "bumps": [
                float(F(edge["delta"])) for edge in family["edges"]
            ],
            "best_games_float": [
                [float(F(value)) for value in game]
                for game in family["games"]
            ],
        }
    if record.get("best_games_float") is None:
        raise ValueError("guided parent must contain floating games")
    return payload, record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attachment", type=int, required=True)
    parser.add_argument("--template-anchor", type=int, required=True)
    parser.add_argument("--permutation", required=True)
    parser.add_argument("--scale", type=float, required=True)
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--envelope-nodes")
    parser.add_argument("--envelope-coalitions")
    parser.add_argument("--envelope-multiplier", type=float, default=4.0)
    parser.add_argument("--envelope-gap", type=float, default=0.0)
    args = parser.parse_args()

    permutation = tuple(
        int(value) for value in args.permutation.split(",")
    )
    if sorted(permutation) != list(range(5)):
        raise ValueError("permutation must contain 0,1,2,3,4")
    _, parent = load_payload(args.parent)
    parent_edges = [
        tuple(int(value) for value in edge) for edge in parent["edges"]
    ]
    edges = list(parent_edges)
    bumps = [float(value) for value in parent["bumps"]]
    game_count = 1 + max(max(edge[0], edge[1]) for edge in edges)
    template = json.loads(args.template.read_text())["family"]
    node_map = {args.template_anchor: args.attachment}
    for node in range(len(template["games"])):
        if node == args.template_anchor:
            continue
        node_map[node] = game_count
        game_count += 1
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
    if args.envelope_nodes is not None:
        nodes = [
            int(value) for value in args.envelope_nodes.split(",")
        ]
        coalitions = [
            int(value)
            for value in args.envelope_coalitions.split(",")
        ]
        if len(nodes) != len(coalitions):
            raise ValueError("envelope node and coalition counts differ")
        row = {}
        for node, coalition in zip(nodes, coalitions, strict=True):
            column = model.game_index(node, coalition)
            row[column] = row.get(column, 0.0) - 1.0
        model.add_ub(
            row,
            -args.envelope_multiplier * (1.0 + args.envelope_gap),
        )
        model.a_ub = model.sparse(model.ub_rows)

    _, iw, ew, metadata = margin_dual(
        parent["best_games_float"], parent_edges
    )
    games = solve_model(model, flow_objective(model, iw, ew, metadata))
    best_margin = float("inf")
    best_games = None
    trace = []
    for _ in range(args.iterations):
        margin, iw, ew, metadata = margin_dual(games, edges)
        trace.append(margin)
        if margin < best_margin:
            best_margin = margin
            best_games = games
        next_games = solve_model(
            model, flow_objective(model, iw, ew, metadata)
        )
        next_margin, _, _, _ = margin_dual(next_games, edges)
        games = next_games
        if next_margin >= margin - 1e-9:
            trace.append(next_margin)
            if next_margin < best_margin:
                best_margin = next_margin
                best_games = games
            break
    if best_games is None:
        raise RuntimeError("guided search produced no family")
    record = {
        "n": 5,
        "edges": [list(edge) for edge in edges],
        "bumps": bumps,
        "best_margin_float": best_margin,
        "best_games_float": best_games,
        "traces": [trace],
        "facet_count": len(model.facets),
        "variables": model.variable_count,
        "inequalities": len(model.ub_rows),
        "equalities": len(model.eq_rows),
    }
    same_grand = all(
        game[31] == best_games[0][31] for game in best_games
    )
    payload = {
        "status": "guided_chain_extension_float",
        "parent": str(args.parent),
        "template": str(args.template),
        "attachment": args.attachment,
        "template_anchor": args.template_anchor,
        "permutation": list(permutation),
        "scale": args.scale,
        "trace": trace,
        "common_core_gap_float": (
            common_core_gap_float(best_games, 5)
            if same_grand
            else None
        ),
        "record": record,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "margin": best_margin,
                "common_core_gap": payload["common_core_gap_float"],
                "trace": trace,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
