#!/usr/bin/env python3
"""Repeat one oriented nonlarge attachment along a designated output port."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from n5_facet_search import permute_mask, search
from n5_transverse_chain_extend import load_record


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--attachment", type=int, required=True)
    parser.add_argument("--template-anchor", type=int, required=True)
    parser.add_argument("--template-output-port", type=int, required=True)
    parser.add_argument("--permutation", required=True)
    parser.add_argument("--scale", type=float, required=True)
    parser.add_argument("--repeat", type=int, required=True)
    parser.add_argument("--starts", type=int, default=6)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--seed", type=int, default=45001)
    args = parser.parse_args()

    permutation = tuple(
        int(value) for value in args.permutation.split(",")
    )
    if sorted(permutation) != list(range(5)):
        raise ValueError("permutation must contain 0,1,2,3,4")
    if args.template_output_port == args.template_anchor:
        raise ValueError("output port must differ from the anchor")

    parent = load_record(args.parent)
    edges = [
        tuple(int(value) for value in edge) for edge in parent["edges"]
    ]
    bumps = [float(value) for value in parent["bumps"]]
    game_count = 1 + max(max(edge[0], edge[1]) for edge in edges)
    attachment = args.attachment
    template = json.loads(args.template.read_text())["family"]
    ports = [attachment]

    for _ in range(args.repeat):
        node_map = {args.template_anchor: attachment}
        for node in range(len(template["games"])):
            if node == args.template_anchor:
                continue
            node_map[node] = game_count
            game_count += 1
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
        attachment = node_map[args.template_output_port]
        ports.append(attachment)

    record = search(
        edges,
        bumps,
        args.starts,
        args.iterations,
        args.seed,
    )
    payload = {
        "status": "repeated_nonlarge_port_float",
        "parent": str(args.parent),
        "template": str(args.template),
        "initial_attachment": args.attachment,
        "template_anchor": args.template_anchor,
        "template_output_port": args.template_output_port,
        "permutation": list(permutation),
        "scale": args.scale,
        "repeat": args.repeat,
        "ports": ports,
        "record": record,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "games": len(record["best_games_float"]),
                "edges": len(record["edges"]),
                "margin": record["best_margin_float"],
                "ports": ports,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
