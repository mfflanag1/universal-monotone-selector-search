#!/usr/bin/env python3
"""Audit legal exact coordinate moves around mixed-premium witnesses."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Iterable


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
sys.path.insert(0, str(PROJECT_SRC))

from n5_facet_search import load_facets


F = Fraction


def dot(left: Iterable[F], right: Iterable[F]) -> F:
    return sum((a * b for a, b in zip(left, right, strict=True)), F(0))


def move_capacity(
    game: tuple[F, ...],
    coalition: int,
    direction: int,
    facets: list[tuple[F, ...]],
    n: int,
) -> F | None:
    bounds: list[F] = []
    for facet in facets:
        coefficient = direction * facet[coalition]
        if coefficient < 0:
            bounds.append(dot(facet, game) / -coefficient)
    for lower in range(1 << n):
        for player in range(n):
            if lower >> player & 1:
                continue
            upper = lower | (1 << player)
            coefficient = direction * (
                int(upper == coalition) - int(lower == coalition)
            )
            if coefficient < 0:
                bounds.append((game[upper] - game[lower]) / -coefficient)
    return min(bounds) if bounds else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("premium_screen", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.premium_screen.read_text())
    n = 5
    grand = (1 << n) - 1
    raw_facets, _ = load_facets()
    facets = [tuple(F(value) for value in facet) for facet in raw_facets]
    audits = []
    exact_archive = payload.get("status") == "n5_three_level_premium_exact_certificate"
    for result in payload["results"]:
        raw_game = result["witness_game"] if exact_archive else result["game"]
        game = tuple(
            F(str(value)).limit_denominator(1_000_000)
            for value in raw_game
        )
        facet_slacks = [dot(facet, game) for facet in facets]
        monotonicity_slacks = [
            game[coalition | (1 << player)] - game[coalition]
            for coalition in range(grand + 1)
            for player in range(n)
            if not coalition >> player & 1
        ]
        if min(facet_slacks) < 0 or min(monotonicity_slacks) < 0:
            raise RuntimeError("rationalized witness left the monotone exact cone")
        upward = []
        downward = []
        for coalition in range(1, grand):
            up = move_capacity(game, coalition, 1, facets, n)
            down = move_capacity(game, coalition, -1, facets, n)
            if up is None or up > 0:
                upward.append(
                    [coalition, "unbounded" if up is None else str(up)]
                )
            if down is None or down > 0:
                downward.append(
                    [coalition, "unbounded" if down is None else str(down)]
                )
        negative_players = [
            player
            for player, value in enumerate(result["divergence"])
            if F(value) < 0
        ]
        positive_players = [
            player
            for player, value in enumerate(result["divergence"])
            if F(value) > 0
        ]
        downward_cover = {
            player: [
                coalition
                for coalition, _ in downward
                if coalition >> player & 1
            ]
            for player in negative_players
        }
        upward_cover = {
            player: [
                coalition
                for coalition, _ in upward
                if coalition >> player & 1
            ]
            for player in positive_players
        }
        audits.append(
            {
                "level_multiplicities": result["level_multiplicities"],
                "premium": str(
                    F(str(result["premium"])).limit_denominator(1_000_000)
                ),
                "game": [str(value) for value in game],
                "minimum_exact_facet_slack": str(min(facet_slacks)),
                "minimum_monotonicity_slack": str(
                    min(monotonicity_slacks)
                ),
                "upward_moves": upward,
                "downward_moves": downward,
                "negative_player_downward_cover": downward_cover,
                "positive_player_upward_cover": upward_cover,
            }
        )
    output = {
        "status": "n5_mixed_premium_witness_bump_audit",
        "audits": audits,
    }
    rendered = json.dumps(output, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
