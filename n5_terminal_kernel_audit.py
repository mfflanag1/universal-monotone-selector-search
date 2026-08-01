#!/usr/bin/env python3
"""Reduce an exact dual terminal-path report to its bipartite cycle kernel."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict, deque
from pathlib import Path
from typing import Any


Vertex = tuple[str, int]
Edge = frozenset[Vertex]


def vertex_record(
    vertex: Vertex,
    terminals: dict[int, dict[str, Any]],
    core_degree: int | None = None,
) -> dict[str, Any]:
    kind, node = vertex
    terminal = terminals[node]
    record = {
        "kind": "source" if kind == "s" else "sink",
        "node": node,
        "coalition": terminal.get("coalition"),
        "scaled_divergence": terminal["scaled_divergence"],
    }
    if core_degree is not None:
        record["core_degree"] = core_degree
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    terminals = {
        int(terminal["node"]): terminal for terminal in report["terminals"]
    }
    adjacency: dict[Vertex, set[Vertex]] = defaultdict(set)
    colors: dict[Edge, set[int]] = defaultdict(set)
    for path in report["terminal_path_decomposition"]:
        source = ("s", int(path["source"]))
        sink = ("q", int(path["sink"]))
        adjacency[source].add(sink)
        adjacency[sink].add(source)
        colors[frozenset((source, sink))].add(int(path["player"]) + 1)
    if {node for _, node in adjacency} != set(terminals):
        raise RuntimeError("path endpoints do not equal the terminal set")

    component_count = 0
    seen: set[Vertex] = set()
    for root in adjacency:
        if root in seen:
            continue
        component_count += 1
        seen.add(root)
        queue = [root]
        while queue:
            vertex = queue.pop()
            for neighbor in adjacency[vertex]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)

    degrees = {vertex: len(neighbors) for vertex, neighbors in adjacency.items()}
    leaves = deque(vertex for vertex, degree in degrees.items() if degree <= 1)
    removed: set[Vertex] = set()
    while leaves:
        vertex = leaves.popleft()
        if vertex in removed:
            continue
        removed.add(vertex)
        for neighbor in adjacency[vertex]:
            if neighbor in removed:
                continue
            degrees[neighbor] -= 1
            if degrees[neighbor] == 1:
                leaves.append(neighbor)
    core = set(adjacency) - removed
    core_degree = {
        vertex: sum(neighbor in core for neighbor in adjacency[vertex])
        for vertex in core
    }
    core_edges = {
        frozenset((vertex, neighbor))
        for vertex in core
        for neighbor in adjacency[vertex]
        if neighbor in core
    }
    branches = {
        vertex for vertex, degree in core_degree.items() if degree >= 3
    }

    used: set[Edge] = set()
    segments = []
    for start in sorted(branches):
        for neighbor in sorted(adjacency[start]):
            edge = frozenset((start, neighbor))
            if neighbor not in core or edge in used:
                continue
            path = [start]
            edge_colors = []
            previous = start
            current = neighbor
            used.add(edge)
            while True:
                path.append(current)
                edge_colors.append(
                    sorted(colors[frozenset((previous, current))])
                )
                if current in branches:
                    break
                choices = [
                    candidate
                    for candidate in adjacency[current]
                    if candidate in core and candidate != previous
                ]
                if len(choices) != 1:
                    raise RuntimeError("two-core contains an unclassified vertex")
                previous, current = current, choices[0]
                used.add(frozenset((previous, current)))
            segments.append(
                {
                    "start": vertex_record(path[0], terminals),
                    "end": vertex_record(path[-1], terminals),
                    "length": len(path) - 1,
                    "vertices": [
                        vertex_record(vertex, terminals) for vertex in path
                    ],
                    "edge_player_sets": edge_colors,
                }
            )
    if branches and used != core_edges:
        raise RuntimeError("not every two-core edge entered a kernel segment")

    vertex_count = len(adjacency)
    edge_count = len(colors)
    result = {
        "status": "exact_terminal_transport_kernel_audit",
        "source": str(args.terminal_report),
        "terminal_count": len(terminals),
        "terminal_role_count": vertex_count,
        "mixed_source_sink_node_count": vertex_count - len(terminals),
        "source_count": sum(vertex[0] == "s" for vertex in adjacency),
        "sink_count": sum(vertex[0] == "q" for vertex in adjacency),
        "terminal_pair_count": edge_count,
        "component_count": component_count,
        "cycle_rank": edge_count - vertex_count + component_count,
        "leaf_count": sum(
            len(neighbors) == 1 for neighbors in adjacency.values()
        ),
        "two_core": {
            "vertex_count": len(core),
            "edge_count": len(core_edges),
            "branch_count": len(branches),
            "branches": [
                vertex_record(vertex, terminals, core_degree[vertex])
                for vertex in sorted(branches)
            ],
        },
        "kernel_segment_count": len(segments),
        "kernel_segments": segments,
        "interpretation": (
            "The cycle rank is the number of independent cycles in the "
            "terminal source-sink incidence graph. Repeated leaf deletion "
            "removes tree attachments; degree-two paths in the remaining "
            "two-core are summarized as kernel segments."
        ),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
