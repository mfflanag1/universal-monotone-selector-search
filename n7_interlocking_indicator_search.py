#!/usr/bin/env python3
"""Globally optimize interlocking indicator terminals with witness exactness."""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_interlocking_terminal_search import (
    connected_incidence,
    random_routing,
    scalar_block_transport_exists,
)
from n5_mixed_terminal_search import Topology
from n7_mixed_terminal_witness_search import WitnessExactCone


def optimize_routing(n, sources, sinks, supports):
    grand = (1 << n) - 1
    arcs = tuple(
        (source, len(sources) + sink, supports[source][sink])
        for source in range(len(sources))
        for sink in range(len(sinks))
        if supports[source][sink]
    )
    topology = Topology(
        2 * len(sources), arcs, tuple(() for _ in range(2 * len(sources)))
    )
    cone = WitnessExactCone(n, topology)
    objective = np.zeros(cone.game_variable_count)
    for source, coalition in enumerate(sources):
        objective[cone.game_column(source, coalition)] += 1.0
    for sink, coalition in enumerate(sinks):
        node = len(sources) + sink
        objective[cone.game_column(node, grand ^ coalition)] += 1.0
        objective[cone.game_column(node, grand)] -= 1.0
    point = cone.maximize(objective)
    if point is None:
        return None
    return float(objective @ point[: cone.game_variable_count]), point


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=7)
    parser.add_argument("--k", type=int, default=3)
    parser.add_argument("--samples", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    grand = (1 << args.n) - 1
    rng = random.Random(args.seed)
    tested = 0
    generated = 0
    best = []
    while tested < args.samples:
        generated += 1
        sources, sinks, supports = random_routing(args.k, args.n, rng)
        if any(coalition in (0, grand) for coalition in (*sources, *sinks)):
            continue
        if not connected_incidence(supports):
            continue
        if scalar_block_transport_exists(sources, sinks, args.n):
            continue
        solved = optimize_routing(args.n, sources, sinks, supports)
        if solved is None:
            continue
        tested += 1
        objective, point = solved
        row = {
            "sources": list(sources),
            "sinks": list(sinks),
            "supports": [list(row) for row in supports],
            "objective": objective,
            "games": [
                [
                    float(value)
                    for value in point[
                        node * (grand + 1) : (node + 1) * (grand + 1)
                    ]
                ]
                for node in range(2 * args.k)
            ],
        }
        best.append(row)
        best.sort(key=lambda candidate: candidate["objective"], reverse=True)
        del best[10:]
        print(json.dumps({"tested": tested, "objective": objective}), flush=True)
        if objective > 1e-8:
            break
    payload = {
        "status": (
            "positive_indicator_terminal_found"
            if best and best[0]["objective"] > 1e-8
            else "no_positive_indicator_terminal_found"
        ),
        "n": args.n,
        "k": args.k,
        "tested": tested,
        "generated": generated,
        "best": best,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key != "best"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
