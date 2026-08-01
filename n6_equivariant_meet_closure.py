#!/usr/bin/env python3
"""Add relative-permutation meets to a compact six-player quotient family."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from n5_equivariant_quotient_lp import quotient_family, solve_quotient
from n5_intrinsic_permutation_union import permute_game
from n5_prune_equivariant_dual_support import components, restrict_family
from n6_dummy_lift_equivariant import add_induced_edges, prune


def load_family(path):
    payload = json.loads(path.read_text())
    family = payload["family"]
    games = [
        tuple(Fraction(value) for value in game) for game in family["games"]
    ]
    edges = {tuple(edge) for edge in family["quotient_edges"]}
    return games, edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--permutations", default="1,0,2,3,4,5")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    permutations = [
        tuple(int(value) for value in text.split(","))
        for text in args.permutations.split(";")
    ]
    for permutation in permutations:
        if sorted(permutation) != list(range(6)):
            parser.error("each permutation must list 0,1,2,3,4,5 exactly once")

    games, edges = load_family(args.input)
    source_count = len(games)
    expanded = list(games)
    expanded_set = set(expanded)
    for permutation in permutations:
        right = [permute_game(game, permutation) for game in games]
        for left_game in games:
            for right_game in right:
                meet = tuple(
                    min(left_value, right_value)
                    for left_value, right_value in zip(
                        left_game, right_game, strict=True
                    )
                )
                if meet not in expanded_set:
                    expanded_set.add(meet)
                    expanded.append(meet)
    print(json.dumps({"phase": "meets", "games": len(expanded)}), flush=True)
    raw_edges = add_induced_edges(expanded)
    qgames, induced_edges, stabilizers = quotient_family(expanded, raw_edges, 6)
    if qgames[:source_count] != games:
        raise RuntimeError("canonical source indices changed")
    qedges = set(induced_edges).union(edges)
    print(
        json.dumps(
            {
                "phase": "quotient",
                "games": len(qgames),
                "edges": len(qedges),
            }
        ),
        flush=True,
    )
    best = None
    for nodes in components(len(qgames), qedges):
        local = restrict_family(qgames, qedges, stabilizers, nodes)
        margin = solve_quotient(*local, 6)
        if best is None or margin < best[0]:
            best = (margin, *local)
        if margin < -1e-8:
            break
    if best is None:
        raise RuntimeError("meet closure has no component")
    margin, games, edges, stabilizers = best
    margin, games, edges, stabilizers, pruning = prune(
        games, edges, stabilizers, args.tolerance
    )
    found = margin < -1e-8
    output = {
        "status": (
            "equivariant_obstruction_found"
            if found
            else "no_equivariant_obstruction_found"
        ),
        "source": str(args.input),
        "permutations": [list(permutation) for permutation in permutations],
        "expanded_games": len(qgames),
        "expanded_edges": len(qedges),
        "margin": margin,
        "pruned_games": len(games),
        "pruned_edges": len(edges),
        "pruning": pruning,
        "family": {
            "games": [list(map(str, game)) for game in games],
            "quotient_edges": [list(edge) for edge in sorted(edges)],
        },
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: value for key, value in output.items() if key != "family"}, indent=2))
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
