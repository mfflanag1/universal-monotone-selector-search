#!/usr/bin/env python3
"""Close selected family nodes through a coordinatewise extremum."""

from __future__ import annotations

import argparse
import itertools
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
GRAND = 31


def parse_ints(text: str) -> list[int]:
    return [int(value) for value in text.split(",")]


def load_parent(path: Path) -> tuple[dict, list[list[F]]]:
    payload = json.loads(path.read_text())
    if payload.get("record") is not None:
        record = payload["record"]
        games = [
            [F(value).limit_denominator(100_000_000) for value in game]
            for game in record["best_games_float"]
        ]
        return record, games
    family = payload["family"]
    games = [[F(value) for value in game] for game in family["games"]]
    record = {
        "edges": [
            [
                int(edge["lower"]),
                int(edge["upper"]),
                int(edge["coalition"]),
            ]
            for edge in family["edges"]
        ],
        "bumps": [float(F(edge["delta"])) for edge in family["edges"]],
        "best_games_float": [
            [float(value) for value in game] for game in games
        ],
    }
    return record, games


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--endpoints", required=True)
    parser.add_argument(
        "--mode", choices=("lower", "upper"), default="lower"
    )
    parser.add_argument(
        "--order",
        choices=("coalition", "reverse", "small-first", "large-first"),
        default="coalition",
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

    parent, games_exact = load_parent(args.parent)
    parent_edges = [
        tuple(int(value) for value in edge) for edge in parent["edges"]
    ]
    edges = list(parent_edges)
    bumps = [float(value) for value in parent["bumps"]]
    endpoints = parse_ints(args.endpoints)
    root = games_exact[0]
    displacements = [
        [
            games_exact[endpoint][coalition] - root[coalition]
            for coalition in range(GRAND + 1)
        ]
        for endpoint in endpoints
    ]
    common = [F(0)] * (GRAND + 1)
    for coalition in range(1, GRAND + 1):
        values = [
            displacement[coalition]
            for displacement in displacements
        ]
        common[coalition] = (
            min(values) if args.mode == "lower" else max(values)
        )

    next_node = len(games_exact)
    common_node = next_node
    next_node += 1
    connector_edges = []
    connector_bumps = []
    arms = []
    for endpoint, displacement in zip(
        endpoints, displacements, strict=True
    ):
        differences = []
        for coalition in range(1, GRAND + 1):
            delta = (
                displacement[coalition] - common[coalition]
                if args.mode == "lower"
                else common[coalition] - displacement[coalition]
            )
            if delta < 0:
                raise ValueError("coordinatewise arm has negative bump")
            if delta:
                differences.append((coalition, delta))
        if args.order == "reverse":
            differences.reverse()
        elif args.order == "small-first":
            differences.sort(key=lambda item: (item[1], item[0]))
        elif args.order == "large-first":
            differences.sort(key=lambda item: (-item[1], item[0]))
        arms.append(
            [[coalition, str(delta)] for coalition, delta in differences]
        )
        if args.closure == "path":
            current = common_node if args.mode == "lower" else endpoint
            for index, (coalition, delta) in enumerate(differences):
                last = index == len(differences) - 1
                target = (
                    endpoint if args.mode == "lower" else common_node
                ) if last else next_node
                if not last:
                    next_node += 1
                connector_edges.append((current, target, coalition))
                connector_bumps.append(float(delta))
                current = target
            continue

        full = (1 << len(differences)) - 1
        nodes = {
            0: common_node if args.mode == "lower" else endpoint,
            full: endpoint if args.mode == "lower" else common_node,
        }
        for size in range(1, len(differences)):
            for indices in itertools.combinations(
                range(len(differences)), size
            ):
                subset = sum(1 << index for index in indices)
                nodes[subset] = next_node
                next_node += 1
        for subset in range(full + 1):
            for index, (coalition, delta) in enumerate(differences):
                bit = 1 << index
                if subset & bit:
                    continue
                connector_edges.append(
                    (nodes[subset], nodes[subset | bit], coalition)
                )
                connector_bumps.append(float(delta))

    edges.extend(connector_edges)
    bumps.extend(connector_bumps)
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
        raise RuntimeError("coordinate closure search produced no family")

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
        "status": "coordinate_closure_float",
        "parent": str(args.parent),
        "endpoints": endpoints,
        "common_node": common_node,
        "mode": args.mode,
        "order": args.order,
        "closure": args.closure,
        "arms": arms,
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
                "common_node": common_node,
                "arm_lengths": [len(arm) for arm in arms],
                "trace": trace,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
