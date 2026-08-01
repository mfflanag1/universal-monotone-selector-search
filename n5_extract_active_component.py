#!/usr/bin/env python3
"""Extract the selector component carrying an archive's active dual flow."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    raw_family = payload["family"]
    games = raw_family["games"]
    edges = raw_family["edges"]
    graph = nx.Graph()
    graph.add_nodes_from(range(len(games)))
    graph.add_edges_from(
        (int(edge["lower"]), int(edge["upper"])) for edge in edges
    )
    active_nodes = {
        int(node)
        for active in payload["exact_margin_dual_certificate"][
            "active_inequalities"
        ]
        for row in [active["row"]]
        if row[0] == "monotonicity"
        for node in (row[1], row[2])
    }
    active_components = [
        component
        for component in nx.connected_components(graph)
        if active_nodes.intersection(component)
    ]
    if len(active_components) != 1:
        raise RuntimeError(
            "active dual flow is not carried by exactly one component"
        )
    source_nodes = sorted(active_components[0])
    remap = {
        source_node: node
        for node, source_node in enumerate(source_nodes)
    }
    component_edges = [
        {
            "lower": remap[int(edge["lower"])],
            "upper": remap[int(edge["upper"])],
            "coalition": int(edge["coalition"]),
        }
        for edge in edges
        if int(edge["lower"]) in remap and int(edge["upper"]) in remap
    ]
    result = {
        "status": "active_selector_component_extracted",
        "n": int(payload["n"]),
        "source": str(args.input),
        "source_nodes": source_nodes,
        "family": {
            "games": [games[node] for node in source_nodes],
            "edges": component_edges,
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "node_count": len(source_nodes),
                "edge_count": len(component_edges),
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
