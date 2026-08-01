#!/usr/bin/env python3
"""Search sparse interlocking protected-flow topologies on six exact games."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from n5_mixed_terminal_search import Topology
from n6_mixed_terminal_search import load_facets, optimize


N = 6


def bipartite_topology(
    rng: random.Random, source_count: int, sink_count: int
) -> Topology | None:
    sink_offset = source_count
    node_count = source_count + sink_count
    arc_players: dict[tuple[int, int], int] = {}
    divergence = [[0] * N for _ in range(node_count)]
    for player in range(N):
        path_count = rng.randint(1, min(3, source_count * sink_count))
        paths = rng.sample(
            [
                (source, sink_offset + sink)
                for source in range(source_count)
                for sink in range(sink_count)
            ],
            path_count,
        )
        for source, sink in paths:
            weight = rng.randint(1, 5)
            arc_players[source, sink] = (
                arc_players.get((source, sink), 0) | (1 << player)
            )
            divergence[source][player] += weight
            divergence[sink][player] -= weight
    if not any(len(set(vector)) >= 3 for vector in divergence):
        return None
    if len(arc_players) < source_count + sink_count - 1:
        return None
    return Topology(
        node_count,
        tuple(
            (source, sink, protected)
            for (source, sink), protected in sorted(arc_players.items())
        ),
        tuple(tuple(vector) for vector in divergence),
    )


def layered_topology(rng: random.Random) -> Topology | None:
    # Two sources, two middle nodes, and two sinks. Each player's flow follows
    # one or two complete paths, so every arc has a literal protected flow.
    node_count = 6
    arc_players: dict[tuple[int, int], int] = {}
    divergence = [[0] * N for _ in range(node_count)]
    for player in range(N):
        path_count = rng.randint(1, 2)
        paths = rng.sample(
            [
                (source, middle, sink)
                for source in range(2)
                for middle in range(2, 4)
                for sink in range(4, 6)
            ],
            path_count,
        )
        for source, middle, sink in paths:
            weight = rng.randint(1, 5)
            for lower, upper in ((source, middle), (middle, sink)):
                arc_players[lower, upper] = (
                    arc_players.get((lower, upper), 0) | (1 << player)
                )
                divergence[lower][player] += weight
                divergence[upper][player] -= weight
    if not any(len(set(vector)) >= 3 for vector in divergence):
        return None
    return Topology(
        node_count,
        tuple(
            (lower, upper, protected)
            for (lower, upper), protected in sorted(arc_players.items())
        ),
        tuple(tuple(vector) for vector in divergence),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--shape", choices=("bipartite22", "bipartite23", "layered222"),
        default="bipartite22"
    )
    parser.add_argument("--samples", type=int, default=20)
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument("--iterations", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facets = load_facets()
    rng = random.Random(args.seed)
    rows = []
    tested = 0
    while tested < args.samples:
        if args.shape == "bipartite22":
            topology = bipartite_topology(rng, 2, 2)
        elif args.shape == "bipartite23":
            topology = bipartite_topology(rng, 2, 3)
        else:
            topology = layered_topology(rng)
        if topology is None:
            continue
        tested += 1
        result = optimize(
            topology, facets, rng, args.starts, args.iterations
        )
        if result is None:
            continue
        row = {
            "arcs": [list(arc) for arc in topology.arcs],
            "divergence": [list(vector) for vector in topology.divergence],
            **result,
        }
        rows.append(row)
        rows.sort(key=lambda item: float(item["objective"]), reverse=True)
        del rows[10:]
        print(
            json.dumps({"tested": tested, "objective": row["objective"]}),
            flush=True,
        )
        if float(row["objective"]) > 1e-8:
            break
    payload = {
        "status": (
            "positive_terminal_obstruction_found"
            if rows and float(rows[0]["objective"]) > 1e-8
            else "no_positive_terminal_obstruction_found"
        ),
        "facet_count": len(facets),
        "tested": tested,
        "configuration": vars(args) | {"output": str(args.output)},
        "best": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
