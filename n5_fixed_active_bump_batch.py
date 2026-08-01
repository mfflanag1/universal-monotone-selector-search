#!/usr/bin/env python3
"""Test direct active-coalition successors of a fixed exact family."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

from family_search import common_core_gap_float
from n5_facet_search import load_facets, margin_dual


F = Fraction
N = 5
GRAND = 31


def parse_ints(text: str) -> list[int]:
    return [int(value) for value in text.split(",")]


def maximum_bump(game: list[F], coalition: int) -> F:
    bounds = []
    for facet in load_facets()[0]:
        coefficient = F(facet[coalition])
        if coefficient >= 0:
            continue
        slack = sum(
            F(facet[item]) * game[item]
            for item in range(1, GRAND + 1)
        )
        bounds.append(slack / -coefficient)
    for superset in range(1, GRAND + 1):
        if coalition & ~superset:
            continue
        if superset == coalition:
            continue
        bounds.append(game[superset] - game[coalition])
    if not bounds:
        raise ValueError("coalition bump is unbounded")
    maximum = min(bounds)
    if maximum <= 0:
        raise ValueError("coalition has no positive exact bump")
    return maximum


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sources", required=True)
    parser.add_argument("--coalitions", required=True)
    parser.add_argument(
        "--factors", default="0.25,0.5,0.75,0.95"
    )
    args = parser.parse_args()

    payload = json.loads(args.parent.read_text())
    family = payload["family"]
    games_exact = [
        [F(value) for value in game] for game in family["games"]
    ]
    edges = [
        (
            int(edge["lower"]),
            int(edge["upper"]),
            int(edge["coalition"]),
        )
        for edge in family["edges"]
    ]
    bumps = [float(F(edge["delta"])) for edge in family["edges"]]
    sources = parse_ints(args.sources)
    coalitions = parse_ints(args.coalitions)
    if len(sources) != len(coalitions):
        raise ValueError("source and coalition counts differ")
    factors = [F(value) for value in args.factors.split(",")]
    maxima = [
        maximum_bump(games_exact[source], coalition)
        for source, coalition in zip(
            sources, coalitions, strict=True
        )
    ]

    best = None
    rows = []
    for factor_tuple in itertools.product(
        factors, repeat=len(sources)
    ):
        candidate_games = [list(game) for game in games_exact]
        candidate_edges = list(edges)
        candidate_bumps = list(bumps)
        chosen_bumps = []
        for source, coalition, maximum, factor in zip(
            sources,
            coalitions,
            maxima,
            factor_tuple,
            strict=True,
        ):
            bump = maximum * factor
            successor = list(games_exact[source])
            successor[coalition] += bump
            target = len(candidate_games)
            candidate_games.append(successor)
            candidate_edges.append((source, target, coalition))
            candidate_bumps.append(float(bump))
            chosen_bumps.append(bump)
        games_float = [
            [float(value) for value in game]
            for game in candidate_games
        ]
        margin, _, _, _ = margin_dual(
            games_float, candidate_edges
        )
        row = {
            "factors": [str(value) for value in factor_tuple],
            "bumps": [str(value) for value in chosen_bumps],
            "margin": margin,
        }
        rows.append(row)
        if best is None or margin < best[0]:
            best = (
                margin,
                factor_tuple,
                chosen_bumps,
                games_float,
                candidate_edges,
                candidate_bumps,
            )
            print(row, flush=True)
    if best is None:
        raise RuntimeError("fixed bump batch produced no family")
    (
        margin,
        factor_tuple,
        chosen_bumps,
        games_float,
        candidate_edges,
        candidate_bumps,
    ) = best
    result = {
        "status": "fixed_active_bump_batch_float",
        "parent": str(args.parent),
        "sources": sources,
        "coalitions": coalitions,
        "maximum_bumps_exact": [str(value) for value in maxima],
        "rows": rows,
        "best_factors": [str(value) for value in factor_tuple],
        "best_bumps": [str(value) for value in chosen_bumps],
        "common_core_gap_float": common_core_gap_float(
            games_float, N
        ),
        "record": {
            "n": N,
            "edges": [list(edge) for edge in candidate_edges],
            "bumps": candidate_bumps,
            "best_margin_float": margin,
            "best_games_float": games_float,
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "maximum_bumps_exact": result[
                    "maximum_bumps_exact"
                ],
                "best_factors": result["best_factors"],
                "best_bumps": result["best_bumps"],
                "margin": margin,
                "common_core_gap": result[
                    "common_core_gap_float"
                ],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
