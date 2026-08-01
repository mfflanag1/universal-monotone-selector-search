#!/usr/bin/env python3
"""Audit common-core gaps separately on selector graph components."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import networkx as nx

from family_search import common_core_gap_exact


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.archive.read_text())
    n = int(payload["n"])
    games = [
        tuple(F(value) for value in game)
        for game in payload["family"]["games"]
    ]
    edges = [
        (
            int(edge["lower"]),
            int(edge["upper"]),
            int(edge["coalition"]),
        )
        for edge in payload["family"]["edges"]
    ]
    graph = nx.Graph()
    graph.add_nodes_from(range(len(games)))
    graph.add_edges_from((lower, upper) for lower, upper, _ in edges)
    active_nodes = {
        int(node)
        for active in payload["exact_margin_dual_certificate"][
            "active_inequalities"
        ]
        for row in [active["row"]]
        if row[0] == "monotonicity"
        for node in (row[1], row[2])
    }

    components = []
    for index, raw_nodes in enumerate(
        sorted(
            nx.connected_components(graph),
            key=lambda nodes: (-len(nodes), min(nodes)),
        )
    ):
        nodes = sorted(raw_nodes)
        component_games = [games[node] for node in nodes]
        same_grand = all(
            game[-1] == component_games[0][-1]
            for game in component_games
        )
        component_edges = sum(
            lower in raw_nodes and upper in raw_nodes
            for lower, upper, _ in edges
        )
        components.append(
            {
                "component": index,
                "node_count": len(nodes),
                "edge_count": component_edges,
                "nodes": nodes,
                "contains_active_dual_flow": bool(
                    active_nodes.intersection(raw_nodes)
                ),
                "common_grand": same_grand,
                "common_core_gap": (
                    str(common_core_gap_exact(component_games, n))
                    if same_grand
                    else None
                ),
            }
        )

    result = {
        "status": "component_common_core_gap_audit",
        "source": str(args.archive),
        "global_common_core_gap": payload.get(
            "common_core_budget_gap_exact"
        ),
        "component_count": len(components),
        "active_dual_component_gaps": [
            component["common_core_gap"]
            for component in components
            if component["contains_active_dual_flow"]
        ],
        "components": components,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
