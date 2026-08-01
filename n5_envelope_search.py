#!/usr/bin/env python3
"""Minimize the family margin while preserving a balanced envelope gap."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from family_search import common_core_gap_float
from n5_facet_search import (
    ExactConeModel,
    flow_objective,
    margin_dual,
    random_objective,
    solve_model,
)
from n5_transverse_chain_extend import load_record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--envelope-nodes", required=True)
    parser.add_argument("--envelope-coalitions", required=True)
    parser.add_argument("--envelope-multiplier", type=float, required=True)
    parser.add_argument("--envelope-gap", type=float, default=0.0)
    parser.add_argument("--starts", type=int, default=8)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42001)
    args = parser.parse_args()

    record = load_record(args.input)
    edges = [
        tuple(int(value) for value in edge) for edge in record["edges"]
    ]
    bumps = [float(value) for value in record["bumps"]]
    nodes = [int(value) for value in args.envelope_nodes.split(",")]
    coalitions = [
        int(value) for value in args.envelope_coalitions.split(",")
    ]
    if len(nodes) != len(coalitions):
        raise ValueError("envelope node and coalition counts differ")

    model = ExactConeModel(edges, bumps)
    row: dict[int, float] = {}
    for node, coalition in zip(nodes, coalitions, strict=True):
        column = model.game_index(node, coalition)
        row[column] = row.get(column, 0.0) - 1.0
    model.add_ub(
        row,
        -args.envelope_multiplier * (1.0 + args.envelope_gap),
    )
    model.a_ub = model.sparse(model.ub_rows)

    rng = random.Random(args.seed)
    best_margin = float("inf")
    best_games = None
    traces = []
    for start in range(args.starts):
        games = solve_model(model, random_objective(model, rng))
        trace = []
        for _ in range(args.iterations):
            margin, iw, ew, metadata = margin_dual(games, edges)
            trace.append(margin)
            if margin < best_margin:
                best_margin = margin
                best_games = games
            objective = flow_objective(model, iw, ew, metadata)
            next_games = solve_model(model, objective)
            next_margin, _, _, _ = margin_dual(next_games, edges)
            games = next_games
            if next_margin >= margin - 1e-9:
                trace.append(next_margin)
                if next_margin < best_margin:
                    best_margin = next_margin
                    best_games = games
                break
        traces.append(trace)
        print(
            json.dumps(
                {
                    "start": start,
                    "trace": trace,
                    "best_margin": best_margin,
                }
            ),
            flush=True,
        )
    if best_games is None:
        raise RuntimeError("search produced no game family")
    payload = {
        "status": "envelope_constrained_float",
        "source": str(args.input),
        "envelope": {
            "nodes": nodes,
            "coalitions": coalitions,
            "multiplier": args.envelope_multiplier,
            "gap": args.envelope_gap,
        },
        "record": {
            "n": 5,
            "edges": [list(edge) for edge in edges],
            "bumps": bumps,
            "starts": args.starts,
            "iterations": args.iterations,
            "seed": args.seed,
            "best_margin_float": best_margin,
            "best_games_float": best_games,
            "traces": traces,
            "facet_count": len(model.facets),
            "variables": model.variable_count,
            "inequalities": len(model.ub_rows),
            "equalities": len(model.eq_rows),
        },
        "common_core_gap_float": common_core_gap_float(best_games, 5),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "margin": best_margin,
                "common_core_gap": payload["common_core_gap_float"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
