#!/usr/bin/env python3
"""Rationalize and recanonicalize an adversarial quotient-family candidate."""

from __future__ import annotations

import argparse
import itertools
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
from n5_equivariant_quotient_lp import bulk_canonical_games, solve_quotient
from n5_facet_search import load_facets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--denominator", type=int, default=100_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    source = Path(payload["source"])
    _source_games, source_edges, _source_stabilizers = load_quotient(source)
    games = [
        tuple(
            Fraction(round(value * args.denominator), args.denominator)
            for value in game
        )
        for game in payload["best_games_float"]
    ]

    facets = load_facets()[0]
    minimum_facet_slack = min(
        sum(
            value * Fraction(coefficient)
            for value, coefficient in zip(game, facet, strict=True)
        )
        for game in games
        for facet in facets
    )
    minimum_monotonicity_slack = min(
        game[coalition | (1 << player)] - game[coalition]
        for game in games
        for coalition in range(32)
        for player in range(5)
        if not coalition >> player & 1
    )
    if minimum_facet_slack < 0 or minimum_monotonicity_slack < 0:
        raise RuntimeError("rationalized candidate left the monotone exact cone")

    canonical_games, maps, stabilizers = bulk_canonical_games(
        games, 5, list(itertools.permutations(range(5)))
    )
    canonical_edges = {
        (
            maps[lower][0],
            maps[lower][1][lower_player],
            maps[upper][0],
            maps[upper][1][upper_player],
        )
        for lower, lower_player, upper, upper_player in source_edges
    }
    margins = {
        method: solve_quotient(
            canonical_games,
            canonical_edges,
            stabilizers,
            5,
            method=method,
        )
        for method in ("highs-ipm", "highs-ds")
    }
    output = {
        "status": "rationalized_topology_candidate",
        "source": str(args.input),
        "denominator": args.denominator,
        "minimum_exact_facet_slack": str(minimum_facet_slack),
        "minimum_monotonicity_slack": str(minimum_monotonicity_slack),
        "margins": margins,
        "family": {
            "games": [list(map(str, game)) for game in canonical_games],
            "quotient_edges": [list(edge) for edge in sorted(canonical_edges)],
        },
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: value for key, value in output.items() if key != "family"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
