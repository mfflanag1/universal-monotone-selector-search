#!/usr/bin/env python3
"""Amplify a rational exact family about its nearest additive Boolean games."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_equivariant_boundary_closure import load_quotient
from n5_equivariant_quotient_lp import solve_quotient
from n5_facet_search import load_facets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--factor", type=Fraction, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    games, edges, stabilizers = load_quotient(args.input)
    scaled_games = []
    for game in games:
        base = tuple(Fraction(round(float(value))) for value in game)
        scaled_games.append(
            tuple(
                base_value + args.factor * (value - base_value)
                for value, base_value in zip(game, base, strict=True)
            )
        )

    facets = load_facets()[0]
    minimum_facet_slack = min(
        sum(
            value * Fraction(coefficient)
            for value, coefficient in zip(game, facet, strict=True)
        )
        for game in scaled_games
        for facet in facets
    )
    minimum_monotonicity_slack = min(
        game[coalition | (1 << player)] - game[coalition]
        for game in scaled_games
        for coalition in range(32)
        for player in range(5)
        if not coalition >> player & 1
    )
    if minimum_facet_slack < 0 or minimum_monotonicity_slack < 0:
        raise RuntimeError("scaled family left the monotone exact cone")

    margins = {
        method: solve_quotient(
            scaled_games, edges, stabilizers, 5, method=method
        )
        for method in ("highs-ipm", "highs-ds")
    }
    output = {
        "status": "scaled_additive_tangent_family",
        "source": str(args.input),
        "factor": str(args.factor),
        "minimum_exact_facet_slack": str(minimum_facet_slack),
        "minimum_monotonicity_slack": str(minimum_monotonicity_slack),
        "margins": margins,
        "family": {
            "games": [list(map(str, game)) for game in scaled_games],
            "quotient_edges": [list(edge) for edge in sorted(edges)],
        },
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: value for key, value in output.items() if key != "family"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
