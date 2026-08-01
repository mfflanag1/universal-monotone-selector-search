#!/usr/bin/env python3
"""List exact n=5 archives with non-box dual tax and their core supports."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


F = Fraction
ROOT = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/results"
)


rows = []
for path in ROOT.glob("*.json"):
    try:
        payload = json.loads(path.read_text())
        tax = F(payload["non_atomic_facet_tax_exact"])
        margin = F(payload["max_min_monotonicity_margin_exact"])
        dual = payload["exact_margin_dual_certificate"]
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        continue
    active_core = [
        (
            int(row["row"][1]),
            int(row["row"][2]),
            str(row["weight"]),
        )
        for row in dual.get("active_inequalities", [])
        if row["row"][0] == "core"
    ]
    sizes = sorted(
        {coalition.bit_count() for _, coalition, _ in active_core}
    )
    if tax <= 0 or sizes != [2, 3]:
        continue
    family = payload.get("family", {})
    games = family.get("games", [])
    edges = family.get("edges", [])
    rational_games = [
        [F(value) for value in game] for game in games
    ]
    active_max_slacks = [
        str(
            max(game[coalition] for game in rational_games)
            - rational_games[node][coalition]
        )
        for node, coalition, _ in active_core
    ]
    rows.append(
        {
            "file": path.name,
            "games": len(games),
            "edges": len(edges),
            "margin": float(margin),
            "tax": float(tax),
            "tax_ratio": float(tax / margin) if margin else None,
            "core": active_core,
            "active_max_slacks": active_max_slacks,
            "all_active_are_maxima": all(
                F(slack) == 0 for slack in active_max_slacks
            ),
            "sizes": sizes,
            "grand_edges": sum(
                int(edge["coalition"]) == 31 for edge in edges
            ),
        }
    )

rows.sort(
    key=lambda row: (
        -row["tax_ratio"] if row["tax_ratio"] is not None else 0,
        row["games"],
    )
)
for row in rows[:100]:
    print(json.dumps(row))
