#!/usr/bin/env python3
"""Exactify the non-sharp K2,3 empty-common-core diagnostic witness."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


F = Fraction
N = 5
GRAND = (1 << N) - 1


def unanimity_game(mask: int) -> list[F]:
    return [F(int(coalition & mask == mask)) for coalition in range(GRAND + 1)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.input.read_text())
    masks = tuple(int(mask) for mask in source["masks"])
    if masks != (4, 15, 8, 4, 23, 8):
        raise AssertionError("unexpected K2,3 masks")
    games = [
        [F(value).limit_denominator(10) for value in game]
        for game in source["games"]
    ]
    unanimity_masks = (12, 12, 4, 12, 8)
    if games != [unanimity_game(mask) for mask in unanimity_masks]:
        raise AssertionError("floating games do not reconstruct as unanimity games")

    arcs = (
        (0, 2, 4),
        (0, 3, 15),
        (0, 4, 8),
        (1, 2, 4),
        (1, 3, 23),
        (1, 4, 8),
    )
    for lower, upper, protected in arcs:
        for coalition in range(GRAND + 1):
            if coalition & protected == protected:
                if games[lower][coalition] > games[upper][coalition]:
                    raise AssertionError("directed protected relation failed")
            elif games[lower][coalition] != games[upper][coalition]:
                raise AssertionError("protected invariance failed")

    divergences = [[0] * N for _ in range(5)]
    for lower, upper, protected in arcs:
        for player in range(N):
            if protected >> player & 1:
                divergences[lower][player] += 1
                divergences[upper][player] -= 1
    lower_expectations = [
        min(
            F(divergences[node][player])
            for player in range(N)
            if unanimity_masks[node] >> player & 1
        )
        for node in range(5)
    ]
    if lower_expectations != [F(2), F(1), F(-2), F(-2), F(-2)]:
        raise AssertionError("unexpected exact lower expectations")

    # Sink 2 fixes player 2's payoff at one; sink 4 fixes player 3's payoff
    # at one. Their joint proper-coalition budget is therefore two against
    # grand worth one, so the exact common-core budget gap is one.
    payload = {
        "status": "n5_k23_empty_common_core_witness_exact",
        "masks": list(masks),
        "unanimity_game_masks": list(unanimity_masks),
        "games": [[str(value) for value in game] for game in games],
        "common_core_budget": "2",
        "grand_worth": "1",
        "common_core_gap_exact": "1",
        "terminal_divergences": divergences,
        "terminal_lower_expectations_exact": [
            str(value) for value in lower_expectations
        ],
        "terminal_objective_exact": str(sum(lower_expectations, F(0))),
        "interpretation": (
            "The terminal cores have empty intersection, but this realization "
            "is safely negative for the protected-flow selector objective."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "common_core_gap_exact": payload["common_core_gap_exact"],
                "terminal_objective_exact": payload["terminal_objective_exact"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
