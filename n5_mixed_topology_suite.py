#!/usr/bin/env python3
"""Search hand-designed terminal topologies with mixed edge flows."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Iterable

from n5_mixed_terminal_search import Topology, optimize_topology


Flow = tuple[int, int, int, int, int]
ArcFlow = tuple[int, int, int, Flow]


def topology(node_count: int, arc_flows: Iterable[ArcFlow]) -> Topology:
    rows = list(arc_flows)
    divergence = [[0] * 5 for _ in range(node_count)]
    for lower, upper, _protected, flow in rows:
        for player, weight in enumerate(flow):
            divergence[lower][player] += weight
            divergence[upper][player] -= weight
    return Topology(
        node_count,
        tuple((lower, upper, protected) for lower, upper, protected, _ in rows),
        tuple(tuple(vector) for vector in divergence),
    )


def cases() -> dict[str, Topology]:
    mu12 = (2, 1, 0, 0, 0)
    mu23 = (0, 2, 1, 0, 0)
    mu34 = (0, 0, 2, 1, 0)
    mu45 = (0, 0, 0, 2, 1)
    mu15 = (2, 0, 0, 0, 1)
    return {
        "fork-disjoint": topology(
            3,
            ((0, 1, 3, mu12), (0, 2, 12, mu34)),
        ),
        "fork-overlap": topology(
            3,
            ((0, 1, 3, mu12), (0, 2, 6, mu23)),
        ),
        "fork-cycle5": topology(
            6,
            (
                (0, 1, 3, mu12),
                (0, 2, 6, mu23),
                (0, 3, 12, mu34),
                (0, 4, 24, mu45),
                (0, 5, 17, mu15),
            ),
        ),
        "cross-disjoint": topology(
            4,
            (
                (0, 2, 3, mu12),
                (0, 3, 12, mu34),
                (1, 2, 12, mu34),
                (1, 3, 3, mu12),
            ),
        ),
        "cross-overlap": topology(
            4,
            (
                (0, 2, 3, mu12),
                (0, 3, 6, mu23),
                (1, 2, 6, mu23),
                (1, 3, 3, mu12),
            ),
        ),
        "ring-transport": topology(
            5,
            (
                (0, 2, 3, mu12),
                (0, 3, 6, mu23),
                (1, 3, 12, mu34),
                (1, 4, 24, mu45),
                (2, 4, 17, mu15),
            ),
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", type=int, default=40)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    rows = []
    for name, candidate in cases().items():
        optimum = optimize_topology(
            candidate, 5, rng, args.starts, args.iterations
        )
        rows.append(
            {
                "name": name,
                "arcs": [list(arc) for arc in candidate.arcs],
                "divergence": [list(vector) for vector in candidate.divergence],
                **({} if optimum is None else optimum),
            }
        )
    rows.sort(key=lambda row: float(row.get("objective", -1e100)), reverse=True)
    payload = {
        "status": "mixed_topology_suite_float",
        "positive_cases": sum(float(row.get("objective", 0.0)) > 1e-8 for row in rows),
        "rows": rows,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
