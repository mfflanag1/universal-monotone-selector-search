#!/usr/bin/env python3
"""Report exact one-coalition bump capacity at selected family nodes."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from n5_fixed_active_bump_batch import maximum_bump


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--coalitions", required=True)
    parser.add_argument("--envelope-maxima-only", action="store_true")
    args = parser.parse_args()

    payload = json.loads(args.archive.read_text())
    games = [
        [F(value) for value in game]
        for game in payload["family"]["games"]
    ]
    rows = []
    for coalition in (
        int(value) for value in args.coalitions.split(",")
    ):
        maximum_value = max(game[coalition] for game in games)
        for node, game in enumerate(games):
            if (
                args.envelope_maxima_only
                and game[coalition] != maximum_value
            ):
                continue
            try:
                capacity = maximum_bump(game, coalition)
            except ValueError:
                capacity = F(0)
            rows.append(
                {
                    "node": node,
                    "coalition": coalition,
                    "envelope_maximum": game[coalition]
                    == maximum_value,
                    "value": str(game[coalition]),
                    "capacity": str(capacity),
                }
            )
    rows.sort(key=lambda row: F(row["capacity"]), reverse=True)
    for row in rows:
        print(json.dumps(row))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
