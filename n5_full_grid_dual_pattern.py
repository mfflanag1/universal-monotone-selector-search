#!/usr/bin/env python3
"""Display active full-grid dual divergences in intrinsic grid coordinates."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import lcm
from pathlib import Path

from exactified_obstruction_paths import exact_closure
from n5_exactified_obstruction_paths import GRAND, null_lift


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.archive.read_text())
    games = [
        tuple(F(value) for value in game)
        for game in payload["family"]["games"]
    ]
    dual = payload["exact_margin_dual_certificate"]
    base4 = tuple(
        map(
            F,
            [0, 0, 0, 0, 0, 3, 3, 1, 0, 3, 3, 1, 3, 1, 5, 6],
        )
    )
    endpoints4 = [list(base4), list(base4)]
    endpoints4[0][7] += 3
    endpoints4[1][14] += 1
    closures = [
        null_lift(exact_closure(game, 4))
        for game in (
            base4,
            tuple(endpoints4[0]),
            tuple(endpoints4[1]),
        )
    ]
    base = closures[0]
    differences = [
        tuple(
            coalition
            for coalition in range(1, GRAND)
            if endpoint[coalition] != base[coalition]
        )
        for endpoint in closures[1:]
    ]

    denominators = [
        F(active["weight"]).denominator
        for active in dual["active_inequalities"]
        if active["row"][0] == "monotonicity"
    ]
    scale = lcm(*denominators)
    divergence = [[F(0)] * 5 for _ in games]
    for active in dual["active_inequalities"]:
        row = active["row"]
        if row[0] != "monotonicity":
            continue
        _, lower, upper, _, player = row
        flow = -F(active["weight"])
        divergence[lower][player] += flow
        divergence[upper][player] -= flow

    records = []
    for node, vector in enumerate(divergence):
        if not any(vector):
            continue
        game = games[node]
        branch = 0
        state: tuple[int, ...] = ()
        for candidate, endpoint in enumerate(closures[1:], start=1):
            candidate_state = tuple(
                int(
                    (game[coalition] - base[coalition])
                    * args.steps
                    / (endpoint[coalition] - base[coalition])
                )
                for coalition in differences[candidate - 1]
            )
            reconstructed = list(base)
            for level, coalition in zip(
                candidate_state,
                differences[candidate - 1],
                strict=True,
            ):
                reconstructed[coalition] += (
                    endpoint[coalition] - base[coalition]
                ) * F(level, args.steps)
            if all(
                reconstructed[coalition] == game[coalition]
                for coalition in range(1, GRAND)
            ):
                branch = candidate
                state = candidate_state
                break
        records.append(
            {
                "node": node,
                "branch": branch,
                "state": state,
                "scaled_divergence": [
                    str(value * scale) for value in vector
                ],
            }
        )
    result = {
        "steps": args.steps,
        "scale": scale,
        "objective_scaled": str(
            -F(dual["objective"]) * scale
        ),
        "active_divergence_nodes": len(records),
        "records": records,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
