#!/usr/bin/env python3
"""Print compact diagnostics for an n=5 floating-point family record."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from family_search import box_family_slack_float, common_core_gap_float
from n5_facet_search import margin_dual


def load_record(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("record") is not None:
        return payload["record"]
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", type=Path, nargs="+")
    args = parser.parse_args()

    for path in args.inputs:
        record = load_record(path)
        games = record["best_games_float"]
        same_grand = all(
            game[31] == games[0][31] for game in games
        )
        edges = [
            tuple(int(value) for value in edge)
            for edge in record["edges"]
        ]
        margin, inequality_weights, _, metadata = margin_dual(games, edges)
        box_margin, _ = box_family_slack_float(games, edges, 5)
        active_core = [
            [list(tag), -weight]
            for weight, tag in zip(
                inequality_weights, metadata, strict=True
            )
            if tag[0] == "core" and abs(weight) > 1e-9
        ]
        print(
            json.dumps(
                {
                    "file": str(path),
                    "games": len(games),
                    "edges": len(edges),
                    "margin": margin,
                    "box_margin": box_margin,
                    "tax": (
                        box_margin - margin
                        if box_margin is not None
                        else None
                    ),
                    "common_core_gap": (
                        common_core_gap_float(games, 5)
                        if same_grand
                        else None
                    ),
                    "active_core": active_core,
                }
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
