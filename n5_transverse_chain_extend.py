#!/usr/bin/env python3
"""Attach one permuted braided block to an existing n=5 float record."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from n5_facet_search import permute_mask, search


F = Fraction


def load_record(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("record") is not None:
        return payload["record"]
    if payload.get("search") is not None:
        return payload["search"]
    if payload.get("family") is not None:
        family = payload["family"]
        normalizer = F(family["games"][0][31])
        return {
            "edges": [
                [
                    int(edge["lower"]),
                    int(edge["upper"]),
                    int(edge["coalition"]),
                ]
                for edge in family["edges"]
            ],
            "bumps": [
                float(F(edge["delta"]) / normalizer)
                for edge in family["edges"]
            ],
        }
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attachment", type=int, required=True)
    parser.add_argument("--template-anchor", type=int, default=38)
    parser.add_argument("--permutation", default="2,1,4,0,3")
    parser.add_argument("--scale", type=float, default=0.00005)
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--seed", type=int, default=15402)
    args = parser.parse_args()

    permutation = tuple(
        int(value) for value in args.permutation.split(",")
    )
    if sorted(permutation) != list(range(5)):
        raise ValueError("permutation must contain 0,1,2,3,4 exactly once")

    parent = load_record(args.parent)
    edges = [tuple(int(value) for value in edge) for edge in parent["edges"]]
    bumps = [float(value) for value in parent["bumps"]]
    game_count = 1 + max(
        max(lower, upper) for lower, upper, _ in edges
    )
    if not 0 <= args.attachment < game_count:
        raise ValueError("attachment is outside the parent game range")

    template = json.loads(args.template.read_text())["family"]
    template_count = len(template["games"])
    if not 0 <= args.template_anchor < template_count:
        raise ValueError("template anchor is outside the template")

    mapping: dict[int, int] = {args.template_anchor: args.attachment}
    next_node = game_count
    for node in range(template_count):
        if node == args.template_anchor:
            continue
        mapping[node] = next_node
        next_node += 1

    for edge in template["edges"]:
        edges.append(
            (
                mapping[int(edge["lower"])],
                mapping[int(edge["upper"])],
                permute_mask(int(edge["coalition"]), permutation),
            )
        )
        bumps.append(args.scale * float(F(edge["delta"])))

    record = search(
        edges,
        bumps,
        args.starts,
        args.iterations,
        args.seed,
    )
    payload = {
        "attachment": args.attachment,
        "template_anchor": args.template_anchor,
        "permutation": list(permutation),
        "scale": args.scale,
        "parent": str(args.parent),
        "template": str(args.template),
        "record": record,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": "saved",
                "games": len(record["best_games_float"]),
                "edges": len(record["edges"]),
                "margin": record["best_margin_float"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
