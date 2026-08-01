#!/usr/bin/env python3
"""Report the active floating margin-dual support of an n=5 family."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from n5_facet_search import margin_dual


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--threshold", type=float, default=1e-9)
    args = parser.parse_args()

    for path in args.inputs:
        payload = json.loads(path.read_text())
        record = payload.get("record", payload)
        if record.get("best_games_float") is None and payload.get("family"):
            family = payload["family"]
            record = {
                "edges": [
                    [
                        int(edge["lower"]),
                        int(edge["upper"]),
                        int(edge["coalition"]),
                    ]
                    for edge in family["edges"]
                ],
                "best_games_float": [
                    [float(Fraction(value)) for value in game]
                    for game in family["games"]
                ],
            }
        edges = [tuple(int(value) for value in edge) for edge in record["edges"]]
        margin, inequality_weights, equality_weights, metadata = margin_dual(
            record["best_games_float"], edges
        )
        active = [
            {"row": list(tag), "weight": -weight}
            for weight, tag in zip(inequality_weights, metadata, strict=True)
            if abs(weight) > args.threshold
        ]
        active_efficiency = [
            {"node": node, "weight": weight}
            for node, weight in enumerate(equality_weights)
            if abs(weight) > args.threshold
        ]
        print(
            json.dumps(
                {
                    "file": str(path),
                    "margin": margin,
                    "active": active,
                    "active_efficiency": active_efficiency,
                }
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
