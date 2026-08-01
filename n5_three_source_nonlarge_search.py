#!/usr/bin/env python3
"""Search cyclically interlocked routings of three Biswas-type sources."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from n5_mixed_terminal_search import Topology, optimize_topology


N = 5


def candidate(rng: random.Random) -> Topology | None:
    ordered_pairs = rng.sample(
        [(first, second) for first in range(N) for second in range(N) if first != second],
        3,
    )
    arc_players: dict[tuple[int, int], int] = {}
    divergence = [[0] * N for _ in range(6)]
    sink_loads = [0, 0, 0]
    for source, (heavy, light) in enumerate(ordered_pairs):
        tokens = (heavy, heavy, light)
        destinations = [rng.randrange(3) for _ in tokens]
        if len(set(destinations)) == 1:
            return None
        for player, sink in zip(tokens, destinations, strict=True):
            upper = 3 + sink
            arc_players[source, upper] = (
                arc_players.get((source, upper), 0) | (1 << player)
            )
            divergence[source][player] += 1
            divergence[upper][player] -= 1
            sink_loads[sink] += 1
    if any(load == 0 for load in sink_loads):
        return None
    adjacency = [set() for _ in range(6)]
    for lower, upper in arc_players:
        adjacency[lower].add(upper)
        adjacency[upper].add(lower)
    if any(len(adjacency[node]) < 2 for node in range(6)):
        return None
    seen = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for neighbor in adjacency[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    if len(seen) != 6:
        return None
    return Topology(
        6,
        tuple(
            (lower, upper, protected)
            for (lower, upper), protected in sorted(arc_players.items())
        ),
        tuple(tuple(row) for row in divergence),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--starts", type=int, default=4)
    parser.add_argument("--iterations", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    tested = 0
    generated = 0
    best = []
    while tested < args.samples:
        generated += 1
        topology = candidate(rng)
        if topology is None:
            continue
        tested += 1
        optimum = optimize_topology(
            topology, N, rng, args.starts, args.iterations
        )
        if optimum is None:
            continue
        row = {
            "arcs": [list(arc) for arc in topology.arcs],
            "divergence": [list(vector) for vector in topology.divergence],
            **optimum,
        }
        best.append(row)
        best.sort(key=lambda item: float(item["objective"]), reverse=True)
        del best[20:]
        if tested % 10 == 0 or float(row["objective"]) > 1e-8:
            print(
                json.dumps(
                    {
                        "tested": tested,
                        "objective": row["objective"],
                        "best": best[0]["objective"],
                    }
                ),
                flush=True,
            )
        if float(row["objective"]) > 1e-8:
            break
    payload = {
        "status": (
            "positive_three_source_topology_found"
            if best and float(best[0]["objective"]) > 1e-8
            else "no_positive_three_source_topology_found"
        ),
        "tested": tested,
        "generated": generated,
        "best": best,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key != "best"}))
    return 1 if payload["status"] == "positive_three_source_topology_found" else 0


if __name__ == "__main__":
    raise SystemExit(main())
