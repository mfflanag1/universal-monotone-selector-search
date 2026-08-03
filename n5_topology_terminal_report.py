#!/usr/bin/env python3
"""Build a terminal report for hand-specified protected block topologies."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from n5_mixed_terminal_search import (
    complete_bipartite_rectangular_topology,
    complete_bipartite_topology,
    weighted_bipartite_six_cycle_topology,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shape", choices=("six-cycle", "k23", "k33"), required=True)
    parser.add_argument("--masks", type=int, nargs="+", required=True)
    parser.add_argument("--weights", type=int, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    expected = {"six-cycle": 6, "k23": 6, "k33": 9}[args.shape]
    if len(args.masks) != expected:
        parser.error(f"{args.shape} requires {expected} masks")
    weights = tuple(args.weights) if args.weights is not None else (1,) * expected
    if len(weights) != expected:
        parser.error(f"{args.shape} requires {expected} weights")
    topology = (
        weighted_bipartite_six_cycle_topology(tuple(args.masks), weights, 5)
        if args.shape == "six-cycle"
        else complete_bipartite_rectangular_topology(
            tuple(args.masks), weights, 2, 3, 5
        )
        if args.shape == "k23"
        else complete_bipartite_topology(tuple(args.masks), weights, 3, 5)
    )
    terminals = [
        {
            "node": node,
            "scaled_divergence": [str(value) for value in divergence],
            "signed_indicator": len(set(divergence)) == 2,
            "coefficient": None,
            "coalition": None,
        }
        for node, divergence in enumerate(topology.divergence)
        if any(divergence)
    ]
    paths = []
    for (source, sink, mask), weight in zip(topology.arcs, weights, strict=True):
        for player in range(5):
            if mask >> player & 1:
                paths.append(
                    {
                        "player": player,
                        "source": source,
                        "sink": sink,
                        "scaled_flow": str(weight),
                        "path": [source, sink],
                    }
                )
    args.output.write_text(
        json.dumps(
            {
                "status": "exact_dual_terminal_report",
                "terminals": terminals,
                "terminal_path_decomposition": paths,
            },
            indent=2,
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
