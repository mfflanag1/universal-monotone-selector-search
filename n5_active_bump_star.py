#!/usr/bin/env python3
"""Extend selected active nodes by direct coalition-bump successors."""

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
    solve_model,
)


F = Fraction


def parse_ints(text: str) -> list[int]:
    return [int(value) for value in text.split(",")]


def parse_floats(text: str) -> list[float]:
    return [float(value) for value in text.split(",")]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sources", required=True)
    parser.add_argument("--coalitions", required=True)
    parser.add_argument("--bumps", required=True)
    parser.add_argument(
        "--direction",
        choices=("outgoing", "incoming"),
        default="outgoing",
    )
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--envelope-nodes", required=True)
    parser.add_argument("--envelope-coalitions", required=True)
    parser.add_argument("--envelope-multiplier", type=float, default=4.0)
    parser.add_argument("--envelope-gap", type=float, default=0.0)
    args = parser.parse_args()

    payload = json.loads(args.parent.read_text())
    parent = payload.get("record", payload)
    if parent.get("best_games_float") is None:
        family = payload["family"]
        normalizer = F(family["games"][0][31])
        parent = {
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
            "best_games_float": [
                [float(F(value) / normalizer) for value in game]
                for game in family["games"]
            ],
        }
    parent_edges = [
        tuple(int(value) for value in edge) for edge in parent["edges"]
    ]
    edges = list(parent_edges)
    bumps = [float(value) for value in parent["bumps"]]
    sources = parse_ints(args.sources)
    coalitions = parse_ints(args.coalitions)
    new_bumps = parse_floats(args.bumps)
    if not (
        len(sources) == len(coalitions) == len(new_bumps)
    ):
        raise ValueError("source, coalition, and bump counts differ")

    next_node = len(parent["best_games_float"])
    added_nodes = []
    for source, coalition, bump in zip(
        sources, coalitions, new_bumps, strict=True
    ):
        if bump <= 0:
            raise ValueError("bumps must be positive")
        edge = (
            (source, next_node, coalition)
            if args.direction == "outgoing"
            else (next_node, source, coalition)
        )
        edges.append(edge)
        bumps.append(bump)
        added_nodes.append(next_node)
        next_node += 1

    model = ExactConeModel(edges, bumps)
    envelope_nodes = parse_ints(args.envelope_nodes)
    envelope_coalitions = parse_ints(args.envelope_coalitions)
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
        raise RuntimeError("active bump search produced no family")

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
    grand_values = [game[31] for game in best_games]
    common_gap = (
        common_core_gap_float(best_games, 5)
        if max(grand_values) - min(grand_values) <= 1e-9
        else None
    )
    result = {
        "status": "active_bump_star_float",
        "parent": str(args.parent),
        "sources": sources,
        "coalitions": coalitions,
        "bumps": new_bumps,
        "direction": args.direction,
        "added_nodes": added_nodes,
        "trace": trace,
        "common_core_gap_float": common_gap,
        "record": record,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "margin": best_margin,
                "common_core_gap": result["common_core_gap_float"],
                "direction": args.direction,
                "added_nodes": added_nodes,
                "trace": trace,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
