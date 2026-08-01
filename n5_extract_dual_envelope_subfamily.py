#!/usr/bin/env python3
"""Extract active dual nodes plus selected common-core envelope nodes."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from family_search import (
    Family,
    box_family_slack_float,
    common_core_gap_float,
    family_slack_float,
    serial_family,
)


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--extra-nodes", required=True)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    raw = payload["family"]
    all_games = [
        tuple(F(value) for value in game) for game in raw["games"]
    ]
    nodes = {
        node
        for active in payload["exact_margin_dual_certificate"][
            "active_inequalities"
        ]
        for row in [active["row"]]
        if row[0] == "monotonicity"
        for node in (int(row[1]), int(row[2]))
    }
    nodes.update(int(value) for value in args.extra_nodes.split(","))
    source_nodes = sorted(nodes)
    games = [all_games[node] for node in source_nodes]
    grand = (1 << int(payload["n"])) - 1
    edges = []
    for left in range(len(games)):
        for right in range(left + 1, len(games)):
            differences = [
                coalition
                for coalition, (left_value, right_value) in enumerate(
                    zip(games[left], games[right], strict=True)
                )
                if left_value != right_value
            ]
            if len(differences) != 1 or differences[0] in (0, grand):
                continue
            coalition = differences[0]
            if games[left][coalition] < games[right][coalition]:
                edges.append((left, right, coalition))
            else:
                edges.append((right, left, coalition))
    family = Family(
        games,
        edges,
        [int(raw["depths"][node]) for node in source_nodes],
        "dual_envelope_induced_subfamily",
    )
    n = int(payload["n"])
    margin, _ = family_slack_float(games, edges, n)
    box_margin, _ = box_family_slack_float(games, edges, n)
    common_gap = common_core_gap_float(games, n)
    result = {
        "status": "dual_envelope_induced_subfamily",
        "n": n,
        "source": str(args.input),
        "source_nodes": source_nodes,
        "node_count": len(games),
        "edge_count": len(edges),
        "common_core_budget_gap_float": common_gap,
        "max_min_monotonicity_margin_float": margin,
        "box_max_min_monotonicity_margin_float": box_margin,
        "non_atomic_facet_tax_float": box_margin - margin,
        "family": serial_family(family, margin, common_gap),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "node_count",
                    "edge_count",
                    "common_core_budget_gap_float",
                    "max_min_monotonicity_margin_float",
                    "box_max_min_monotonicity_margin_float",
                    "non_atomic_facet_tax_float",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
