#!/usr/bin/env python3
"""Search additive exact generators whose meets tighten an equivariant family."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from fractions import Fraction
from pathlib import Path

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_equivariant_boundary_closure import best_component, load_quotient, prune
from n5_equivariant_quotient_lp import quotient_family, solve_quotient
from n5_intrinsic_permutation_union import permute_game


def additive_game(point):
    return tuple(
        sum(point[player] for player in range(5) if coalition >> player & 1)
        for coalition in range(32)
    )


def compositions(total, parts=5):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for suffix in compositions(total - first, parts - 1):
            yield (first, *suffix)


def solve_extension(source_games, source_edges, generator, permutations, tolerance):
    expanded = list(source_games)
    expanded_set = set(expanded)
    for permutation in permutations:
        right = permute_game(generator, permutation)
        for game in source_games:
            meet = tuple(
                min(left_value, right_value)
                for left_value, right_value in zip(game, right, strict=True)
            )
            if meet not in expanded_set:
                expanded_set.add(meet)
                expanded.append(meet)
    raw_edges = set()
    add_induced_edges_fast(expanded, raw_edges)
    qgames, induced_edges, stabilizers = quotient_family(
        expanded, raw_edges, 5
    )
    if qgames[: len(source_games)] != source_games:
        raise RuntimeError("canonical source indices changed")
    qedges = set(induced_edges).union(source_edges)
    best = best_component(qgames, qedges, stabilizers)
    if best is None:
        raise RuntimeError("additive meet family has no comparison component")
    margin, games, edges, stabilizers = best
    margin, games, edges, stabilizers, pruning = prune(
        games, edges, stabilizers, tolerance
    )
    return margin, games, edges, len(qgames), len(qedges), pruning


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--total", type=int, default=20)
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260801)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source_games, source_edges, source_stabilizers = load_quotient(args.input)
    source_margin = solve_quotient(
        source_games, source_edges, source_stabilizers, 5
    )
    grand = source_games[0][-1]
    pool = list(compositions(args.total))
    rng = random.Random(args.seed)
    rng.shuffle(pool)
    seeds = [
        (args.total,) * 0 + tuple(
            args.total // 5 for _ in range(5)
        )
    ] if args.total % 5 == 0 else []
    if args.total % 9 == 0:
        scale = args.total // 9
        seeds.extend(
            [
                tuple(scale * value for value in (3, 1, 1, 2, 2)),
                tuple(scale * value for value in (2, 3, 2, 1, 1)),
            ]
        )
    points = list(dict.fromkeys(seeds + pool))[: args.trials]
    permutations = list(itertools.permutations(range(5)))
    best = (source_margin, source_games, source_edges, None, None)
    rows = []
    for trial, counts in enumerate(points):
        point = tuple(grand * Fraction(value, args.total) for value in counts)
        generator = additive_game(point)
        margin, games, edges, qgame_count, qedge_count, pruning = solve_extension(
            source_games,
            source_edges,
            generator,
            permutations,
            args.tolerance,
        )
        row = {
            "trial": trial,
            "counts": list(counts),
            "margin": margin,
            "expanded_games": qgame_count,
            "expanded_edges": qedge_count,
            "pruned_games": len(games),
            "pruned_edges": len(edges),
        }
        rows.append(row)
        if margin < best[0] - 1e-12:
            best = (margin, games, edges, counts, pruning)
            print(json.dumps({"phase": "best", **row}), flush=True)
        elif trial % 10 == 0:
            print(json.dumps(row), flush=True)
        checkpoint = {
            "status": "additive_meet_search_checkpoint",
            "source": str(args.input),
            "source_margin": source_margin,
            "best_margin": best[0],
            "best_counts": best[3],
            "rows": rows,
            "pruning": best[4],
            "family": {
                "games": [list(map(str, game)) for game in best[1]],
                "quotient_edges": [list(edge) for edge in sorted(best[2])],
            },
        }
        args.output.write_text(json.dumps(checkpoint, indent=2) + "\n")
        if margin < -1e-8:
            break
    found = best[0] < -1e-8
    checkpoint["status"] = (
        "equivariant_obstruction_found"
        if found
        else "no_equivariant_obstruction_found"
    )
    args.output.write_text(json.dumps(checkpoint, indent=2) + "\n")
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
