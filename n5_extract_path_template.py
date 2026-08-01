#!/usr/bin/env python3
"""Extract an induced shortest-path exact-family template."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from n5_mixed_product_path import load_archive, shortest_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sources", required=True)
    parser.add_argument("--targets", required=True)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    games, edges = load_archive(args.input)
    path = shortest_path(
        len(games),
        edges,
        [int(value) for value in args.sources.split(",")],
        {int(value) for value in args.targets.split(",")},
    )
    path_set = set(path)
    node_map = {
        old_node: new_node for new_node, old_node in enumerate(path)
    }
    family = {
        "shape": "induced_shortest_path",
        "games": [
            [str(value) for value in games[old_node]]
            for old_node in path
        ],
        "edges": [
            {
                "lower": node_map[lower],
                "upper": node_map[upper],
                "coalition": coalition,
                "delta": str(delta),
            }
            for lower, upper, coalition, delta in edges
            if lower in path_set and upper in path_set
        ],
    }
    result = {
        "status": "exact_path_template",
        "n": int(payload["n"]),
        "source": str(args.input),
        "source_path": path,
        "source_to_template_node": {
            str(old_node): new_node
            for old_node, new_node in node_map.items()
        },
        "node_count": len(path),
        "edge_count": len(family["edges"]),
        "family": family,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "node_count": result["node_count"],
                "edge_count": result["edge_count"],
                "source_start": path[0],
                "source_end": path[-1],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
