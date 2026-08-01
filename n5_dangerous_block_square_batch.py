#!/usr/bin/env python3
"""Search the 22 serious no-common-player five-player block squares.

These are the permutation/row/column representatives with four distinct
incomparable size-two/three protected blocks, full player support, and the
sharp 1/2 local premium at all four K2,2 terminals.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from n5_mixed_terminal_search import block_square_topology, optimize_topology


ALL_MAX_PREMIUM_ORBITS = (
    (1, 7, 7, 2),
    (1, 7, 7, 10),
    (1, 7, 7, 14),
    (1, 7, 11, 2),
    (1, 7, 11, 12),
    (1, 7, 11, 14),
    (1, 7, 11, 18),
    (1, 7, 15, 2),
    (1, 7, 15, 10),
    (1, 7, 25, 10),
    (1, 7, 27, 2),
    (1, 7, 27, 10),
    (1, 15, 7, 10),
    (1, 15, 15, 2),
    (1, 15, 15, 6),
    (1, 15, 19, 6),
    (1, 15, 23, 2),
    (1, 15, 23, 6),
    (3, 5, 5, 6),
    (3, 5, 5, 12),
    (3, 5, 5, 14),
    (3, 5, 5, 28),
    (3, 5, 6, 11),
    (3, 5, 6, 12),
    (3, 5, 6, 15),
    (3, 5, 6, 28),
    (3, 5, 9, 14),
    (3, 5, 9, 28),
    (3, 5, 10, 12),
    (3, 5, 10, 15),
    (3, 5, 10, 19),
    (3, 5, 10, 22),
    (3, 5, 10, 25),
    (3, 5, 10, 28),
    (3, 5, 13, 14),
    (3, 5, 13, 28),
    (3, 5, 14, 22),
    (3, 5, 14, 28),
    (3, 5, 25, 28),
    (3, 5, 26, 28),
    (3, 13, 5, 14),
    (3, 13, 5, 28),
    (3, 13, 6, 21),
    (3, 13, 6, 28),
    (3, 13, 13, 6),
    (3, 13, 13, 14),
    (3, 13, 13, 20),
    (3, 13, 13, 28),
    (3, 13, 14, 3),
    (3, 13, 14, 7),
    (3, 13, 14, 20),
    (3, 13, 14, 28),
    (3, 13, 15, 6),
    (3, 13, 17, 28),
    (3, 13, 18, 28),
    (3, 13, 21, 6),
    (3, 13, 21, 28),
    (3, 13, 22, 3),
    (3, 13, 22, 7),
    (3, 13, 22, 10),
    (3, 13, 22, 14),
    (3, 13, 22, 21),
    (3, 13, 22, 24),
    (3, 13, 22, 28),
    (3, 13, 23, 6),
    (3, 13, 23, 20),
    (3, 15, 15, 12),
    (7, 11, 13, 14),
)


SERIOUS_CANDIDATES = (
    (3, 5, 6, 28),
    (3, 5, 9, 28),
    (3, 5, 10, 22),
    (3, 5, 10, 25),
    (3, 5, 10, 28),
    (3, 5, 14, 22),
    (3, 5, 14, 28),
    (3, 5, 25, 28),
    (3, 5, 26, 28),
    (3, 13, 6, 21),
    (3, 13, 6, 28),
    (3, 13, 14, 20),
    (3, 13, 14, 28),
    (3, 13, 17, 28),
    (3, 13, 18, 28),
    (3, 13, 21, 6),
    (3, 13, 21, 28),
    (3, 13, 22, 10),
    (3, 13, 22, 14),
    (3, 13, 22, 21),
    (3, 13, 22, 24),
    (3, 13, 22, 28),
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", type=int, default=8)
    parser.add_argument("--iterations", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260801)
    parser.add_argument("--all-orbits", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    candidates = (
        ALL_MAX_PREMIUM_ORBITS if args.all_orbits else SERIOUS_CANDIDATES
    )
    rows = []
    positive = 0
    for index, masks in enumerate(candidates):
        topology = block_square_topology(masks, 5)
        result = optimize_topology(
            topology,
            5,
            random.Random(args.seed + index),
            args.starts,
            args.iterations,
        )
        if result is None:
            raise RuntimeError(f"optimization failed for masks {masks}")
        row = {
            "masks": list(masks),
            "arcs": [list(arc) for arc in topology.arcs],
            "divergence": [list(vector) for vector in topology.divergence],
            **result,
        }
        positive += float(result["objective"]) > 1e-8
        rows.append(row)
        print(
            json.dumps(
                {
                    "index": index,
                    "masks": masks,
                    "objective": result["objective"],
                }
            ),
            flush=True,
        )
    rows.sort(key=lambda row: float(row["objective"]), reverse=True)
    payload = {
        "status": "dangerous_block_square_coordinate_ascent_float",
        "candidate_count": len(candidates),
        "all_orbits": args.all_orbits,
        "positive_count": positive,
        "starts": args.starts,
        "iterations": args.iterations,
        "results": rows,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key != "results"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
