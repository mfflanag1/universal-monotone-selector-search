#!/usr/bin/env python3
"""Rank exact archives by genuinely mixed active dual divergence."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("--limit", type=int, default=100)
    args = parser.parse_args()

    rows = []
    for path in args.results.glob("*.json"):
        try:
            payload = json.loads(path.read_text())
            n = int(payload["n"])
            margin = F(payload["max_min_monotonicity_margin_exact"])
            tax = F(payload["non_atomic_facet_tax_exact"])
            dual = payload["exact_margin_dual_certificate"]
            game_count = len(payload["family"]["games"])
        except (KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue
        divergence = [[F(0)] * n for _ in range(game_count)]
        for active in dual.get("active_inequalities", []):
            row = active["row"]
            if row[0] != "monotonicity":
                continue
            _, lower, upper, _, player = row
            flow = -F(active["weight"])
            divergence[int(lower)][int(player)] += flow
            divergence[int(upper)][int(player)] -= flow
        mixed = []
        indicator = []
        for node, vector in enumerate(divergence):
            nonzero = {value for value in vector if value}
            if not nonzero:
                continue
            record = {
                "node": node,
                "divergence": [str(value) for value in vector],
            }
            if len(nonzero) == 1:
                indicator.append(record)
            else:
                mixed.append(record)
        if not mixed:
            continue
        rows.append(
            {
                "file": path.name,
                "games": game_count,
                "margin": str(margin),
                "margin_float": float(margin),
                "tax": str(tax),
                "tax_ratio": float(tax / margin) if margin else None,
                "common_gap": payload.get(
                    "common_core_budget_gap_exact"
                ),
                "mixed_nodes": len(mixed),
                "indicator_nodes": len(indicator),
                "mixed": mixed,
            }
        )
    rows.sort(
        key=lambda row: (
            0 if row["common_gap"] not in (None, "0") else 1,
            -row["tax_ratio"] if row["tax_ratio"] is not None else 0,
            row["margin_float"],
        )
    )
    for row in rows[: args.limit]:
        print(json.dumps(row))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
