#!/usr/bin/env python3
"""Close two bridge endpoints through a common lower connector."""

from __future__ import annotations

import argparse
import json
import itertools
from pathlib import Path

from n5_facet_search import (
    ExactConeModel,
    flow_objective,
    margin_dual,
    solve_model,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--endpoint-a", type=int, required=True)
    parser.add_argument("--endpoint-b", type=int, required=True)
    parser.add_argument("--a-coalitions", required=True)
    parser.add_argument("--b-coalitions", required=True)
    parser.add_argument("--bump", type=float, required=True)
    parser.add_argument(
        "--mode", choices=("lower", "upper"), default="lower"
    )
    parser.add_argument(
        "--closure", choices=("path", "cube"), default="path"
    )
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--envelope-nodes", required=True)
    parser.add_argument("--envelope-coalitions", required=True)
    parser.add_argument("--envelope-multiplier", type=float, default=4.0)
    parser.add_argument("--envelope-gap", type=float, default=0.0)
    args = parser.parse_args()

    payload = json.loads(args.parent.read_text())
    parent = payload.get("record", payload)
    parent_edges = [
        tuple(int(value) for value in edge) for edge in parent["edges"]
    ]
    edges = list(parent_edges)
    bumps = [float(value) for value in parent["bumps"]]
    game_count = 1 + max(max(edge[0], edge[1]) for edge in edges)
    lower = game_count
    a_coalitions = [
        int(value) for value in args.a_coalitions.split(",")
    ]
    b_coalitions = [
        int(value) for value in args.b_coalitions.split(",")
    ]
    if not a_coalitions or not b_coalitions:
        raise ValueError("connector arms cannot be empty")
    connector_edges = []
    next_node = game_count + 1
    for endpoint, coalitions in (
        (args.endpoint_a, a_coalitions),
        (args.endpoint_b, b_coalitions),
    ):
        if args.closure == "path":
            current = lower if args.mode == "lower" else endpoint
            for index, coalition in enumerate(coalitions):
                last = index == len(coalitions) - 1
                target = (
                    endpoint if args.mode == "lower" else lower
                ) if last else next_node
                if not last:
                    next_node += 1
                connector_edges.append((current, target, coalition))
                current = target
            continue

        full = (1 << len(coalitions)) - 1
        nodes = {
            0: lower if args.mode == "lower" else endpoint,
            full: endpoint if args.mode == "lower" else lower,
        }
        for size in range(1, len(coalitions)):
            for indices in itertools.combinations(
                range(len(coalitions)), size
            ):
                subset = sum(1 << index for index in indices)
                nodes[subset] = next_node
                next_node += 1
        for subset in range(full + 1):
            for index, coalition in enumerate(coalitions):
                bit = 1 << index
                if subset & bit:
                    continue
                connector_edges.append(
                    (nodes[subset], nodes[subset | bit], coalition)
                )
    edges.extend(connector_edges)
    bumps.extend([args.bump] * len(connector_edges))

    model = ExactConeModel(edges, bumps)
    envelope_nodes = [
        int(value) for value in args.envelope_nodes.split(",")
    ]
    envelope_coalitions = [
        int(value)
        for value in args.envelope_coalitions.split(",")
    ]
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
        raise RuntimeError("connector search produced no family")
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
        "status": "theta_connector_float",
        "parent": str(args.parent),
        "lower": lower,
        "connector_nodes": list(range(game_count, next_node)),
        "endpoint_a": args.endpoint_a,
        "endpoint_b": args.endpoint_b,
        "a_coalitions": a_coalitions,
        "b_coalitions": b_coalitions,
        "bump": args.bump,
        "mode": args.mode,
        "closure": args.closure,
        "trace": trace,
        "record": record,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "margin": best_margin,
                "trace": trace,
                "connector_nodes": result["connector_nodes"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
