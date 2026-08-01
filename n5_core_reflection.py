#!/usr/bin/env python3
"""Reflect exact cores and reconstruct every induced comparison."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_equivariant_boundary_closure import best_component, load_quotient, prune
from n5_equivariant_quotient_lp import quotient_family
from n5_intrinsic_permutation_union import permute_game


def reflect(game, cap):
    grand = game[31]
    reflected_grand = 5 * cap - grand
    if reflected_grand <= 0:
        raise ValueError("reflection cap must exceed one fifth of the grand worth")
    return tuple(
        (
            coalition.bit_count() * cap
            - grand
            + game[31 ^ coalition]
        )
        * grand
        / reflected_grand
        for coalition in range(32)
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--cap", type=Fraction, default=Fraction(1))
    parser.add_argument("--keep-full", action="store_true")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source_games, _source_edges, _source_stabilizers = load_quotient(args.input)
    reflected = [reflect(game, args.cap) for game in source_games]
    permutations = list(itertools.permutations(range(5)))
    orbit_games = list(
        dict.fromkeys(
            permute_game(game, permutation)
            for game in reflected
            for permutation in permutations
        )
    )
    raw_edges = set()
    add_induced_edges_fast(orbit_games, raw_edges)
    games, edges, stabilizers = quotient_family(orbit_games, raw_edges, 5)
    full_games = list(games)
    full_edges = set(edges)
    best = best_component(games, edges, stabilizers)
    if best is None:
        raise RuntimeError("reflected family has no comparison component")
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
        "cap": str(args.cap),
        "orbit_games": len(orbit_games),
        "raw_edges": len(raw_edges),
        "quotient_games": len(games),
        "quotient_edges": len(edges),
        "margin": margin,
        "pruning": pruning,
        "family": {
            "games": [
                list(map(str, game))
                for game in (full_games if args.keep_full else games)
            ],
            "quotient_edges": [
                list(edge)
                for edge in sorted(full_edges if args.keep_full else edges)
            ],
        },
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: value for key, value in output.items() if key != "family"}, indent=2))
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
