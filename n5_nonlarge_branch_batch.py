#!/usr/bin/env python3
"""Screen simultaneous nonlarge-star branches at one parent game."""

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


def unique_orientations(
    template: dict, target_coalition: int
) -> list[tuple[int, ...]]:
    pair_directions = {
        int(edge["coalition"])
        for edge in template["edges"]
        if int(edge["coalition"]).bit_count() == 2
    }
    representatives: dict[tuple[int, ...], tuple[int, ...]] = {}
    for permutation in itertools.permutations(range(5)):
        signature = tuple(
            sorted(
                permute_mask(int(edge["coalition"]), permutation)
                for edge in template["edges"]
            )
        )
        if target_coalition not in {
            permute_mask(coalition, permutation)
            for coalition in pair_directions
        }:
            continue
        representatives.setdefault(signature, permutation)
    return list(representatives.values())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attachment", type=int, required=True)
    parser.add_argument("--target-coalition", type=int, required=True)
    parser.add_argument("--template-anchor", type=int, default=0)
    parser.add_argument("--branch-count", type=int, default=2)
    parser.add_argument("--scale", type=float, default=0.00001)
    parser.add_argument("--starts", type=int, default=4)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--seed", type=int, default=31001)
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
    orientations = unique_orientations(template, args.target_coalition)
    if not 1 <= args.branch_count <= len(orientations):
        raise ValueError("branch count is outside the orientation range")

    rows = []
    best_margin = None
    best_gap = None
    best_empty = None
    combinations = list(
        itertools.combinations(orientations, args.branch_count)
    )
    for index, combination in enumerate(combinations):
        edges = list(parent_edges)
        bumps = list(parent_bumps)
        next_node = game_count
        for permutation in combination:
            node_map = {args.template_anchor: args.attachment}
            for node in range(len(template["games"])):
                if node == args.template_anchor:
                    continue
                node_map[node] = next_node
                next_node += 1
            for edge in template["edges"]:
                edges.append(
                    (
                        node_map[int(edge["lower"])],
                        node_map[int(edge["upper"])],
                        permute_mask(
                            int(edge["coalition"]), permutation
                        ),
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
            "permutations": [list(value) for value in combination],
            "margin": record["best_margin_float"],
            "common_core_gap": gap,
        }
        rows.append(row)
        candidate = (
            record["best_margin_float"],
            gap,
            combination,
            record,
        )
        if best_margin is None or candidate[0] < best_margin[0]:
            best_margin = candidate
            print(json.dumps({"best_margin": row}), flush=True)
        if best_gap is None or (
            gap is not None and gap > best_gap[1]
        ):
            best_gap = candidate
            print(json.dumps({"best_gap": row}), flush=True)
        if gap is not None and gap > 0 and (
            best_empty is None or candidate[0] < best_empty[0]
        ):
            best_empty = candidate
            print(json.dumps({"best_empty_core": row}), flush=True)

    def serialize(candidate: tuple | None) -> dict | None:
        if candidate is None:
            return None
        return {
            "margin": candidate[0],
            "common_core_gap": candidate[1],
            "permutations": [
                list(value) for value in candidate[2]
            ],
            "record": candidate[3],
        }

    payload = {
        "status": "nonlarge_branch_batch_float",
        "parent": str(args.parent),
        "template": str(args.template),
        "attachment": args.attachment,
        "target_coalition": args.target_coalition,
        "branch_count": args.branch_count,
        "scale": args.scale,
        "orientation_count": len(orientations),
        "rows": rows,
        "best_margin": serialize(best_margin),
        "best_gap": serialize(best_gap),
        "best_empty_core": serialize(best_empty),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "orientation_count": len(orientations),
                "combination_count": len(combinations),
                "best_margin": payload["best_margin"]["margin"],
                "best_margin_gap": payload["best_margin"][
                    "common_core_gap"
                ],
                "best_gap": payload["best_gap"]["common_core_gap"],
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
