#!/usr/bin/env python3
"""List exact n=5 empty-common-core archives and margin-dual supports."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


F = Fraction
ROOTS = [
    Path(
        "/Users/maxf/projects/economics-research/game-theory/"
        "exact-game-monotone-selection/results"
    ),
    Path("/private/tmp"),
]


rows = []
seen = set()
for root in ROOTS:
    for path in root.glob("*.json"):
        try:
            payload = json.loads(path.read_text())
            gap = F(payload["common_core_budget_gap_exact"])
            margin = F(payload["max_min_monotonicity_margin_exact"])
            dual = payload["exact_margin_dual_certificate"]
        except (KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue
        if gap <= 0 or path.name in seen:
            continue
        seen.add(path.name)
        active_core = [
            (int(row["row"][1]), int(row["row"][2]))
            for row in dual.get("active_inequalities", [])
            if row["row"][0] == "core"
        ]
        sizes = sorted(
            {coalition.bit_count() for _, coalition in active_core}
        )
        family = payload.get("family", {})
        rows.append(
            {
                "file": str(path),
                "games": len(family.get("games", [])),
                "edges": len(family.get("edges", [])),
                "gap": float(gap),
                "margin": float(margin),
                "core": active_core,
                "sizes": sizes,
            }
        )

rows.sort(key=lambda row: (row["margin"], -row["gap"]))
for row in rows:
    print(json.dumps(row))
