#!/usr/bin/env python3
"""Glue both endpoints of an exact path template to existing family nodes."""

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


def load_record(path: Path) -> dict:
    payload = json.loads(path.read_text())
    record = payload.get("record", payload)
    if record.get("best_games_float") is None:
        raise ValueError("parent must contain floating games")
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=int, required=True)
    parser.add_argument("--target", type=int, required=True)
    parser.add_argument("--template-source", type=int, required=True)
    parser.add_argument("--template-target", type=int, required=True)
    parser.add_argument("--permutation", required=True)
    parser.add_argument("--scale", type=float, required=True)
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--envelope-nodes", required=True)
    parser.add_argument("--envelope-coalitions", required=True)
    parser.add_argument("--envelope-multiplier", type=float, default=4.0)
    parser.add_argument("--envelope-gap", type=float, default=0.0)
    args = parser.parse_args()

    permutation = tuple(int(value) for value in args.permutation.split(","))
    if sorted(permutation) != list(range(5)):
        raise ValueError("permutation must contain 0,1,2,3,4")

    parent = load_record(args.parent)
    parent_edges = [
        tuple(int(value) for value in edge) for edge in parent["edges"]
    ]
    edges = list(parent_edges)
    bumps = [float(value) for value in parent["bumps"]]
    game_count = len(parent["best_games_float"])

    template = json.loads(args.template.read_text())["family"]
    node_map = {
        args.template_source: args.source,
        args.template_target: args.target,
    }
    for node in range(len(template["games"])):
        if node in node_map:
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
    envelope_nodes = [
        int(value) for value in args.envelope_nodes.split(",")
    ]
    envelope_coalitions = [
        int(value) for value in args.envelope_coalitions.split(",")
    ]
    if len(envelope_nodes) != len(envelope_coalitions):
        raise ValueError("envelope node and coalition counts differ")
    row = {}
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
        raise RuntimeError("two-end bridge search produced no family")

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
    result = {
        "status": "two_end_bridge_float",
        "parent": str(args.parent),
        "template": str(args.template),
        "source": args.source,
        "target": args.target,
        "template_source": args.template_source,
        "template_target": args.template_target,
        "permutation": list(permutation),
        "scale": args.scale,
        "trace": trace,
        "common_core_gap_float": common_core_gap_float(best_games, 5),
        "record": record,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "games": len(best_games),
                "edges": len(edges),
                "margin": best_margin,
                "common_core_gap": result["common_core_gap_float"],
                "trace": trace,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
