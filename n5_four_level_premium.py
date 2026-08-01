#!/usr/bin/env python3
"""Screen sharp premiums for five-player four-level unit-step divergences."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from n5_mixed_divergence_premium import maximize_branch_premium, sparse_matrix
from n5_mixed_terminal_branch_search import extreme_representations
from n5_facet_search import load_facets


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    n = 5
    grand = (1 << n) - 1
    facets, _ = load_facets()
    rows: list[dict[int, float]] = []
    rhs: list[float] = []
    for facet in facets:
        rows.append(
            {
                coalition: -float(coefficient)
                for coalition, coefficient in enumerate(facet)
                if coefficient
            }
        )
        rhs.append(0.0)
    for coalition in range(grand + 1):
        for player in range(n):
            if coalition >> player & 1:
                continue
            rows.append({coalition: 1.0, coalition | (1 << player): -1.0})
            rhs.append(0.0)
    equalities = np.zeros((2, grand + 1))
    equalities[0, 0] = 1.0
    equalities[1, grand] = 1.0
    inequality_matrix = sparse_matrix(rows, grand + 1)
    inequality_rhs = np.asarray(rhs)
    results = []
    for counts in ((1, 1, 1, 2), (1, 1, 2, 1), (1, 2, 1, 1), (2, 1, 1, 1)):
        divergence = tuple(
            value
            for level, count in enumerate(counts)
            for value in (F(level),) * count
        )
        representations = extreme_representations(list(divergence), n)
        best = None
        for branch, representation in enumerate(representations):
            solved = maximize_branch_premium(
                divergence,
                representation,
                n,
                inequality_matrix,
                inequality_rhs,
                equalities,
                np.asarray([0.0, 1.0]),
            )
            if solved is None:
                continue
            premium, game = solved
            if best is None or premium > best["premium"]:
                best = {
                    "premium": premium,
                    "branch": branch,
                    "game": [float(value) for value in game],
                }
        if best is None:
            raise RuntimeError("no feasible premium branch")
        row = {
            "level_multiplicities": list(counts),
            "divergence": [str(value) for value in divergence],
            "representation_count": len(representations),
            **best,
        }
        results.append(row)
        print(json.dumps({key: row[key] for key in row if key != "game"}), flush=True)
    args.output.write_text(
        json.dumps(
            {"status": "n5_four_level_premium_float_screen", "results": results},
            indent=2,
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
