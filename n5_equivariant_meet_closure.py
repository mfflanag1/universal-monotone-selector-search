#!/usr/bin/env python3
"""Add relative-permutation meets to a compact equivariant exact family."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_equivariant_boundary_closure import best_component, load_quotient, prune
from n5_equivariant_quotient_lp import quotient_family, solve_quotient
from n5_intrinsic_permutation_union import permute_game
from n5_prune_equivariant_dual_support import components, restrict_family


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--right-input", type=Path)
    parser.add_argument(
        "--permutations",
        default="0,1,2,3,4",
        help="semicolon-separated relative permutations",
    )
    parser.add_argument("--all-permutations", action="store_true")
    parser.add_argument("--source-component-only", action="store_true")
    parser.add_argument("--top-components", type=int, default=0)
    parser.add_argument("--right-alpha", type=Fraction, default=Fraction(1))
    parser.add_argument("--right-point", default="1/5,1/5,1/5,1/5,1/5")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    permutations = (
        list(itertools.permutations(range(5)))
        if args.all_permutations
        else [
            tuple(int(value) for value in text.split(","))
            for text in args.permutations.split(";")
        ]
    )
    for permutation in permutations:
        if sorted(permutation) != list(range(5)):
            parser.error("each permutation must list 0,1,2,3,4 exactly once")

    left_games, edges, _stabilizers = load_quotient(args.input)
    if args.right_input is None:
        right_source_games = left_games
        games = list(left_games)
    else:
        right_source_games, right_edges, _right_stabilizers = load_quotient(
            args.right_input
        )
        games = list(left_games)
        game_index = {game: index for index, game in enumerate(games)}
        right_remap = []
        for game in right_source_games:
            target = game_index.get(game)
            if target is None:
                target = len(games)
                game_index[game] = target
                games.append(game)
            right_remap.append(target)
        edges = set(edges)
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
    point = tuple(Fraction(value) for value in args.right_point.split(","))
    if len(point) != 5 or sum(point) != 1:
        parser.error("--right-point must contain five fractions summing to one")
    if not 0 <= args.right_alpha <= 1:
        parser.error("--right-alpha must lie in [0,1]")
    grand = games[0][-1]
    if any(game[-1] != grand for game in right_source_games):
        raise ValueError("left and right games must have the same grand worth")
    additive = tuple(
        grand
        * sum(point[player] for player in range(5) if coalition >> player & 1)
        for coalition in range(32)
    )
    right_games = [
        tuple(
            args.right_alpha * value
            + (1 - args.right_alpha) * additive_value
            for value, additive_value in zip(game, additive, strict=True)
        )
        for game in right_source_games
    ]
    expanded = list(games)
    expanded_set = set(expanded)
    for permutation in permutations:
        right = [permute_game(game, permutation) for game in right_games]
        for left_game in left_games:
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

    raw_edges = set()
    add_induced_edges_fast(expanded, raw_edges)
    qgames, induced_qedges, stabilizers = quotient_family(
        expanded, raw_edges, 5
    )
    if qgames[: len(games)] != games:
        raise RuntimeError("canonical source indices changed")
    qedges = set(induced_qedges).union(edges)
    ranked_components = []
    if args.source_component_only or args.top_components:
        best = None
        for nodes in components(len(qgames), qedges):
            if args.source_component_only and not any(
                node < source_count for node in nodes
            ):
                continue
            local = restrict_family(qgames, qedges, stabilizers, nodes)
            local_margin = solve_quotient(*local, 5)
            ranked_components.append((local_margin, local))
            if best is None or local_margin < best[0]:
                best = (local_margin, *local)
    else:
        best = best_component(qgames, qedges, stabilizers)
    if best is None:
        raise RuntimeError("meet expansion has no component")
    margin, games, edges, stabilizers = best
    margin, games, edges, stabilizers, pruning = prune(
        games, edges, stabilizers, args.tolerance
    )
    found = margin < -1e-8
    top_components = []
    for local_margin, local in sorted(
        ranked_components, key=lambda item: item[0]
    )[: args.top_components]:
        (
            compact_margin,
            compact_games,
            compact_edges,
            _compact_stabilizers,
            compact_pruning,
        ) = prune(*local, args.tolerance)
        top_components.append(
            {
                "margin": compact_margin,
                "original_margin": local_margin,
                "pruning": compact_pruning,
                "family": {
                    "games": [
                        list(map(str, game)) for game in compact_games
                    ],
                    "quotient_edges": [
                        list(edge) for edge in sorted(compact_edges)
                    ],
                },
            }
        )
    output = {
        "status": (
            "equivariant_obstruction_found"
            if found
            else "no_equivariant_obstruction_found"
        ),
        "source": str(args.input),
        "right_source": str(args.right_input or args.input),
        "permutations": [list(permutation) for permutation in permutations],
        "right_alpha": str(args.right_alpha),
        "right_point": [str(value) for value in point],
        "expanded_games": len(qgames),
        "expanded_edges": len(qedges),
        "margin": margin,
        "pruned_games": len(games),
        "pruned_edges": len(edges),
        "pruning": pruning,
        "top_components": top_components,
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
