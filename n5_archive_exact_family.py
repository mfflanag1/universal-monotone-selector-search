#!/usr/bin/env python3
"""Create an exact archive from a rational family payload."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from n5_exact_archive import exact_archive


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shape", required=True)
    parser.add_argument("--max-denominator", type=int, default=1_000_000)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    family = payload["family"]
    games = [
        tuple(F(value) for value in game) for game in family["games"]
    ]
    edges = [
        (
            int(edge["lower"]),
            int(edge["upper"]),
            int(edge["coalition"]),
        )
        for edge in family["edges"]
    ]
    archive = exact_archive(
        games,
        edges,
        args.shape,
        args.max_denominator,
    )
    archive["source"] = str(args.input)
    archive["source_nodes"] = payload.get("source_nodes")
    args.output.write_text(json.dumps(archive, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: archive[key]
                for key in (
                    "status",
                    "node_count",
                    "edge_count",
                    "common_core_budget_gap_exact",
                    "max_min_monotonicity_margin_exact",
                    "box_max_min_monotonicity_margin_exact",
                    "non_atomic_facet_tax_exact",
                )
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
