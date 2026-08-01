#!/usr/bin/env python3
"""Reconstruct a glued float family from one rational root and exact offsets."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from n5_exactify_constrained_topology import (
    glued_archive_offsets,
    parse_ints,
)
from n5_facet_search import exact_archive


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--gadget", type=Path, required=True)
    parser.add_argument("--attachment", type=int, required=True)
    parser.add_argument("--gadget-anchor", type=int, required=True)
    parser.add_argument("--gadget-scale", required=True)
    parser.add_argument("--gadget-permutation", default="0,1,2,3,4")
    parser.add_argument("--float-record", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-denominator", type=int, default=100_000_000)
    args = parser.parse_args()

    base = json.loads(args.base.read_text())
    gadget = json.loads(args.gadget.read_text())
    permutation = parse_ints(args.gadget_permutation)
    offsets, edges = glued_archive_offsets(
        base["family"],
        gadget["family"],
        args.attachment,
        args.gadget_anchor,
        F(args.gadget_scale),
        permutation,
    )
    payload = json.loads(args.float_record.read_text())
    record = payload.get("record", payload)
    record_edges = [
        tuple(int(value) for value in edge)
        for edge in record["edges"]
    ]
    if record_edges != edges:
        raise ValueError("float and exact glued topologies differ")
    root = [
        F(value).limit_denominator(args.max_denominator)
        for value in record["best_games_float"][0]
    ]
    exact_games = [
        [
            str(root[coalition] + offset[coalition])
            for coalition in range(len(root))
        ]
        for offset in offsets
    ]
    record["best_games_exact"] = exact_games
    archive = exact_archive(
        record, "n5_glued_root_reconstruction", args.max_denominator
    )
    archive["glued_root_reconstruction"] = {
        "base": str(args.base),
        "gadget": str(args.gadget),
        "attachment": args.attachment,
        "gadget_anchor": args.gadget_anchor,
        "gadget_scale": str(F(args.gadget_scale)),
        "gadget_permutation": list(permutation),
        "float_record": str(args.float_record),
        "maximum_root_denominator": max(
            value.denominator for value in root
        ),
    }
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
