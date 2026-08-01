#!/usr/bin/env python3
"""Add pointwise maxima that remain inside the five-player exact cone."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_equivariant_boundary_closure import best_component, load_quotient, prune
from n5_equivariant_quotient_lp import quotient_family
from n5_facet_search import load_facets
from n5_intrinsic_permutation_union import permute_game


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--right-input", type=Path)
    parser.add_argument(
        "--permutations",
        default="0,1,2,3,4;1,0,2,3,4;1,0,3,2,4;1,2,0,3,4;"
        "1,2,0,4,3;1,2,3,0,4;1,2,3,4,0",
    )
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--exact-root", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    permutations = [
        tuple(int(value) for value in text.split(","))
        for text in args.permutations.split(";")
    ]
    left, left_edges, _left_stabilizers = load_quotient(args.input)
    right, right_edges, _right_stabilizers = load_quotient(
        args.right_input or args.input
    )
    games = list(left)
    game_index = {game: index for index, game in enumerate(games)}
    right_remap = []
    for game in right:
        target = game_index.get(game)
        if target is None:
            target = len(games)
            game_index[game] = target
            games.append(game)
        right_remap.append(target)
    edges = set(left_edges)
    edges.update(
        (
            right_remap[lower],
            lower_player,
            right_remap[upper],
            upper_player,
        )
        for lower, lower_player, upper, upper_player in right_edges
    )
    source_count = len(games)

    candidates = set()
    for permutation in permutations:
        permuted = [permute_game(game, permutation) for game in right]
        candidates.update(
            tuple(
                max(left_value, right_value)
                for left_value, right_value in zip(
                    left_game, right_game, strict=True
                )
            )
            for left_game in left
            for right_game in permuted
        )
    candidates.difference_update(game_index)
    facets = load_facets()[0]
    facet_matrix = np.asarray(facets, dtype=float)
    zero_columns = [
        index for index, facet in enumerate(facets) if not facet[31]
    ]
    rooted_columns = [
        index for index, facet in enumerate(facets) if facet[31]
    ]
    candidate_list = list(candidates)
    accepted = []
    for start in range(0, len(candidate_list), 2000):
        batch = candidate_list[start : start + 2000]
        values = np.asarray(
            [[float(value) for value in game] for game in batch]
        )
        products = values @ facet_matrix.T
        if args.exact_root:
            valid = np.min(products[:, zero_columns], axis=1) >= -1e-10
        else:
            valid = np.min(products, axis=1) >= -1e-10
        for game, product_row, is_valid in zip(
            batch, products, valid, strict=True
        ):
            if not is_valid:
                continue
            if args.exact_root:
                close_zero_columns = [
                    index
                    for index in zero_columns
                    if product_row[index] < 1e-7
                ]
                if any(
                    sum(
                        value * coefficient
                        for value, coefficient in zip(
                            game, facet, strict=True
                        )
                    )
                    < 0
                    for facet in (facets[index] for index in close_zero_columns)
                ):
                    continue
                root = max(game[-1], max(game[:-1]))
                float_root = max(
                    float(root),
                    max(
                        float(game[-1])
                        - product_row[index] / facets[index][31]
                        for index in rooted_columns
                    ),
                )
                active_root_columns = [
                    index
                    for index in rooted_columns
                    if float_root
                    - (
                        float(game[-1])
                        - product_row[index] / facets[index][31]
                    )
                    < 1e-7
                ]
                for index in active_root_columns:
                    facet = facets[index]
                    proper = sum(
                        game[coalition] * facet[coalition]
                        for coalition in range(31)
                    )
                    root = max(root, -proper / facet[31])
                rooted = list(game)
                rooted[31] = root
                accepted.append(tuple(rooted))
            elif all(
                sum(
                    value * coefficient
                    for value, coefficient in zip(game, facet, strict=True)
                )
                >= 0
                for facet in facets
            ):
                accepted.append(game)
    games.extend(accepted)
    print(
        json.dumps(
            {
                "phase": "maxima",
                "candidates": len(candidate_list),
                "accepted_exact": len(accepted),
                "games": len(games),
            }
        ),
        flush=True,
    )
    raw_edges = set()
    add_induced_edges_fast(games, raw_edges)
    qgames, induced_edges, stabilizers = quotient_family(games, raw_edges, 5)
    if qgames[:source_count] != games[:source_count]:
        raise RuntimeError("canonical source indices changed")
    qedges = set(induced_edges).union(edges)
    best = best_component(qgames, qedges, stabilizers)
    if best is None:
        raise RuntimeError("maximum closure has no component")
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
        "right_source": str(args.right_input or args.input),
        "permutations": [list(permutation) for permutation in permutations],
        "candidate_maxima": len(candidate_list),
        "accepted_exact_maxima": len(accepted),
        "exact_root": args.exact_root,
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
