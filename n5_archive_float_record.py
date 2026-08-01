#!/usr/bin/env python3
"""Rationally reconstruct and archive an n=5 floating topology record."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from n5_facet_search import exact_archive


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shape", required=True)
    parser.add_argument("--max-denominator", type=int, default=1_000_000)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    record = payload.get("record", payload)
    archive = exact_archive(
        record,
        args.shape,
        args.max_denominator,
    )
    archive["float_source"] = str(args.input)
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
