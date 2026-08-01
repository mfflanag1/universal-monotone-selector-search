#!/usr/bin/env python3
"""Screen nonlarge-star orientations at a selected recursive attachment port."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

from family_search import common_core_gap_float
from n5_facet_search import permute_mask, search
from n5_transverse_chain_extend import load_record


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attachment", type=int, required=True)
    parser.add_argument("--target-coalition", type=int, required=True)
    parser.add_argument("--template-anchor", type=int, default=0)
    parser.add_argument("--scale", type=float, default=0.00001)
    parser.add_argument("--starts", type=int, default=4)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--seed", type=int, default=30001)
    args = parser.parse_args()

    parent = load_record(args.parent)
    parent_edges = [
        tuple(int(value) for value in edge) for edge in parent["edges"]
    ]
    parent_bumps = [float(value) for value in parent["bumps"]]
    game_count = 1 + max(
        max(lower, upper) for lower, upper, _ in parent_edges
    )
    template = json.loads(args.template.read_text())["family"]
    pair_directions = {
        int(edge["coalition"])
        for edge in template["edges"]
        if int(edge["coalition"]).bit_count() == 2
    }
    permutations = [
        permutation
        for permutation in itertools.permutations(range(5))
        if args.target_coalition
        in {
            permute_mask(coalition, permutation)
            for coalition in pair_directions
        }
    ]

    rows = []
    best_margin = None
    best_gap = None
    for index, permutation in enumerate(permutations):
        node_map = {args.template_anchor: args.attachment}
        next_node = game_count
        for node in range(len(template["games"])):
            if node == args.template_anchor:
                continue
            node_map[node] = next_node
            next_node += 1
        edges = list(parent_edges)
        bumps = list(parent_bumps)
        for edge in template["edges"]:
            edges.append(
                (
                    node_map[int(edge["lower"])],
                    node_map[int(edge["upper"])],
                    permute_mask(int(edge["coalition"]), permutation),
                )
            )
            bumps.append(args.scale * float(F(edge["delta"])))
        record = search(
            edges,
            bumps,
            args.starts,
            args.iterations,
            args.seed + index,
            quiet=True,
        )
        gap = common_core_gap_float(record["best_games_float"], 5)
        row = {
            "permutation": list(permutation),
            "margin": record["best_margin_float"],
            "common_core_gap": gap,
        }
        rows.append(row)
        if best_margin is None or (
            record["best_margin_float"]
            < best_margin[0]
        ):
            best_margin = (
                record["best_margin_float"],
                gap,
                permutation,
                record,
            )
            print(json.dumps({"best_margin": row}), flush=True)
        if gap is not None and gap > 0 and (
            best_gap is None
            or record["best_margin_float"] < best_gap[0]
        ):
            best_gap = (
                record["best_margin_float"],
                gap,
                permutation,
                record,
            )
            print(json.dumps({"best_empty_core": row}), flush=True)

    payload = {
        "status": "nonlarge_orientation_batch_float",
        "parent": str(args.parent),
        "template": str(args.template),
        "attachment": args.attachment,
        "target_coalition": args.target_coalition,
        "scale": args.scale,
        "rows": rows,
        "best_margin": {
            "margin": best_margin[0],
            "common_core_gap": best_margin[1],
            "permutation": list(best_margin[2]),
            "record": best_margin[3],
        },
        "best_empty_core": (
            {
                "margin": best_gap[0],
                "common_core_gap": best_gap[1],
                "permutation": list(best_gap[2]),
                "record": best_gap[3],
            }
            if best_gap is not None
            else None
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "orientation_count": len(permutations),
                "best_margin": payload["best_margin"]["margin"],
                "best_margin_gap": payload["best_margin"][
                    "common_core_gap"
                ],
                "best_empty_core_margin": (
                    payload["best_empty_core"]["margin"]
                    if payload["best_empty_core"] is not None
                    else None
                ),
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
