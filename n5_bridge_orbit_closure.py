#!/usr/bin/env python3
"""Attach bridge orientations and close their endpoints coordinatewise."""

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
    permute_mask,
    solve_model,
)


F = Fraction
GRAND = 31


def parse_permutations(text: str) -> list[tuple[int, ...]]:
    permutations = [
        tuple(int(value) for value in item.split(","))
        for item in text.split(";")
    ]
    if not permutations:
        raise ValueError("at least one permutation is required")
    if any(sorted(permutation) != list(range(5)) for permutation in permutations):
        raise ValueError("each permutation must contain 0,1,2,3,4")
    return permutations


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attachment", type=int, required=True)
    parser.add_argument("--template-anchor", type=int, required=True)
    parser.add_argument("--template-endpoint", type=int, required=True)
    parser.add_argument("--permutations", required=True)
    parser.add_argument("--scale", type=F, required=True)
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

    parent_payload = json.loads(args.parent.read_text())
    parent = parent_payload.get("record", parent_payload)
    parent_edges = [
        tuple(int(value) for value in edge) for edge in parent["edges"]
    ]
    edges = list(parent_edges)
    bumps = [float(value) for value in parent["bumps"]]
    game_count = len(parent["best_games_float"])

    template = json.loads(args.template.read_text())["family"]
    template_games = [
        [F(value) for value in game] for game in template["games"]
    ]
    permutations = parse_permutations(args.permutations)
    endpoints = []
    endpoint_displacements = []
    for permutation in permutations:
        node_map = {args.template_anchor: args.attachment}
        for node in range(len(template_games)):
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
            bumps.append(float(args.scale * F(edge["delta"])))
        endpoints.append(node_map[args.template_endpoint])
        displacement = [F(0)] * (GRAND + 1)
        for coalition in range(1, GRAND + 1):
            mapped = permute_mask(coalition, permutation)
            displacement[mapped] = args.scale * (
                template_games[args.template_endpoint][coalition]
                - template_games[args.template_anchor][coalition]
            )
        endpoint_displacements.append(displacement)

    common = [F(0)] * (GRAND + 1)
    for coalition in range(1, GRAND + 1):
        values = [
            displacement[coalition]
            for displacement in endpoint_displacements
        ]
        common[coalition] = (
            min(values) if args.mode == "lower" else max(values)
        )

    common_node = game_count
    game_count += 1
    connector_edges = []
    connector_bumps = []
    arms = []
    for endpoint, displacement in zip(
        endpoints, endpoint_displacements, strict=True
    ):
        differences = []
        for coalition in range(1, GRAND + 1):
            difference = (
                displacement[coalition] - common[coalition]
                if args.mode == "lower"
                else common[coalition] - displacement[coalition]
            )
            if difference < 0:
                raise ValueError("coordinatewise closure has a negative arm")
            if difference:
                differences.append((coalition, difference))
        arms.append(
            [[coalition, str(delta)] for coalition, delta in differences]
        )

        if args.closure == "path":
            current = common_node if args.mode == "lower" else endpoint
            for index, (coalition, delta) in enumerate(differences):
                last = index == len(differences) - 1
                target = (
                    endpoint if args.mode == "lower" else common_node
                ) if last else game_count
                if not last:
                    game_count += 1
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
                nodes[subset] = game_count
                game_count += 1
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

    envelope_nodes = [
        int(value) for value in args.envelope_nodes.split(",")
    ]
    envelope_coalitions = [
        int(value) for value in args.envelope_coalitions.split(",")
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
        raise RuntimeError("orbit closure search produced no family")

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
        "status": "bridge_orbit_closure_float",
        "parent": str(args.parent),
        "template": str(args.template),
        "attachment": args.attachment,
        "template_anchor": args.template_anchor,
        "template_endpoint": args.template_endpoint,
        "permutations": [list(permutation) for permutation in permutations],
        "scale": str(args.scale),
        "mode": args.mode,
        "closure": args.closure,
        "endpoints": endpoints,
        "common_node": common_node,
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
                "endpoints": endpoints,
                "common_node": common_node,
                "arm_lengths": [len(arm) for arm in arms],
                "trace": trace,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
