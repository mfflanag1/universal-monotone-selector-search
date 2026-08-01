#!/usr/bin/env python3
"""Add every comparison induced by the full permutation closure."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_equivariant_boundary_closure import best_component, load_quotient, prune
from n5_equivariant_quotient_lp import quotient_family
from n5_intrinsic_permutation_union import permute_game


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    games, edges, _stabilizers = load_quotient(args.input)
    permutations = list(itertools.permutations(range(5)))
    orbit_games = list(
        dict.fromkeys(
            permute_game(game, permutation)
            for game in games
            for permutation in permutations
        )
    )
    raw_edges = set()
    add_induced_edges_fast(orbit_games, raw_edges)
    qgames, induced_edges, stabilizers = quotient_family(
        orbit_games, raw_edges, 5
    )
    if qgames[: len(games)] != games:
        raise RuntimeError("canonical source indices changed")
    qedges = set(induced_edges).union(edges)
    best = best_component(qgames, qedges, stabilizers)
    if best is None:
        raise RuntimeError("completed family has no comparison component")
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
        "orbit_games": len(orbit_games),
        "raw_edges": len(raw_edges),
        "quotient_games": len(qgames),
        "quotient_edges": len(qedges),
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
