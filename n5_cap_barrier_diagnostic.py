#!/usr/bin/env python3
"""Diagnose the exact-facet cut blocking an intrinsic coordinate path."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

from exactified_obstruction_paths import exact_closure
from n5_exactified_obstruction_paths import (
    GRAND,
    exact_root,
    null_lift,
)
from n5_facet_search import load_facets


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--grand-cap", type=F, required=True)
    parser.add_argument("--branch", type=int, choices=(1, 2), default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base4 = tuple(
        map(
            F,
            [0, 0, 0, 0, 0, 3, 3, 1, 0, 3, 3, 1, 3, 1, 5, 6],
        )
    )
    endpoints4 = [list(base4), list(base4)]
    endpoints4[0][7] += 3
    endpoints4[1][14] += 1
    base = null_lift(exact_closure(base4, 4))
    endpoint = null_lift(
        exact_closure(tuple(endpoints4[args.branch - 1]), 4)
    )
    differences = tuple(
        coalition
        for coalition in range(1, GRAND)
        if endpoint[coalition] != base[coalition]
    )
    facets, facet_types = load_facets()

    def game_for_state(state: tuple[int, ...]) -> tuple[F, ...]:
        game = list(base)
        for level, coalition in zip(
            state, differences, strict=True
        ):
            game[coalition] = (
                base[coalition]
                + (endpoint[coalition] - base[coalition])
                * F(level, args.steps)
            )
        return tuple(game)

    origin = (0,) * len(differences)
    destination = (args.steps,) * len(differences)
    frontier = {origin}
    dead_record = None
    for depth in range(1, args.steps * len(differences) + 1):
        candidates = {
            tuple(
                level + (index == direction)
                for index, level in enumerate(state)
            )
            for state in frontier
            for direction in range(len(differences))
            if state[direction] < args.steps
        }
        roots = {
            state: exact_root(game_for_state(state), facets)
            for state in candidates
        }
        feasible = {
            state
            for state, root in roots.items()
            if root is not None and root <= args.grand_cap
        }
        if not feasible:
            finite_roots = {
                state: root
                for state, root in roots.items()
                if root is not None
            }
            minimum_root = min(finite_roots.values())
            minimum_states = sorted(
                state
                for state, root in finite_roots.items()
                if root == minimum_root
            )
            active_types = Counter()
            active_facets = Counter()
            for state in minimum_states:
                game = game_for_state(state)
                for index, facet in enumerate(facets):
                    grand_coefficient = facet[GRAND]
                    if grand_coefficient <= 0:
                        continue
                    constant = sum(
                        facet[coalition] * game[coalition]
                        for coalition in range(1, GRAND)
                    )
                    threshold = -constant / grand_coefficient
                    if threshold == minimum_root:
                        active_types[facet_types[index]] += 1
                        active_facets[index] += 1
            dead_record = {
                "dead_depth": depth,
                "previous_frontier_count": len(frontier),
                "candidate_count": len(candidates),
                "finite_root_count": len(finite_roots),
                "minimum_blocking_root": str(minimum_root),
                "minimum_blocking_states": [
                    list(state) for state in minimum_states
                ],
                "active_facet_type_counts": {
                    str(key): value
                    for key, value in sorted(active_types.items())
                },
                "active_facet_indices": sorted(active_facets),
            }
            break
        frontier = feasible
        if destination in frontier:
            break

    result = {
        "status": (
            "blocked_below_cap"
            if dead_record is not None
            else "destination_reached"
        ),
        "steps": args.steps,
        "grand_cap": str(args.grand_cap),
        "branch": args.branch,
        "directions": list(differences),
        "destination_reached": destination in frontier,
        "barrier": dead_record,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
