#!/usr/bin/env python3
"""Lift a five-player equivariant family and close it under S6."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from n5_equivariant_boundary_closure import load_quotient
from n5_equivariant_quotient_lp import solve_quotient
from n5_intrinsic_permutation_union import permute_game
from n5_prune_equivariant_dual_support import components, restrict_family


def dummy_lift(game):
    return tuple(game[coalition & 31] for coalition in range(64))


def veto_lift(game):
    return tuple(
        game[coalition & 31] if coalition >> 5 & 1 else game[0]
        for coalition in range(64)
    )


def add_induced_edges(games):
    edges = set()
    encoded = [
        tuple((value.numerator, value.denominator) for value in game)
        for game in games
    ]
    for coalition in range(1, 63):
        groups = {}
        for index, game in enumerate(encoded):
            signature = game[:coalition] + game[coalition + 1 :]
            groups.setdefault(signature, []).append(index)
        for indices in groups.values():
            if len(indices) < 2:
                continue
            ordered = sorted(indices, key=lambda index: games[index][coalition])
            for position, lower in enumerate(ordered):
                for upper in ordered[position + 1 :]:
                    edges.add((lower, upper, coalition))
    return edges


def dual_games_n6(result, game_count, edges, tolerance):
    core_rows = 62 * game_count
    marginals = np.asarray(result.ineqlin.marginals)
    active = {
        row // 62
        for row in np.flatnonzero(np.abs(marginals[:core_rows]) > tolerance)
    }
    ordered_edges = sorted(edges)
    for edge_index in np.flatnonzero(np.abs(marginals[core_rows:]) > tolerance):
        lower, _lower_player, upper, _upper_player = ordered_edges[int(edge_index)]
        active.update((lower, upper))
    return active


def prune(games, edges, stabilizers, tolerance):
    rows = []
    while True:
        margin, result = solve_quotient(
            games, edges, stabilizers, 6, return_result=True
        )
        active = dual_games_n6(result, len(games), edges, tolerance)
        rows.append(
            {
                "games": len(games),
                "edges": len(edges),
                "dual_games": len(active),
                "margin": margin,
            }
        )
        if not active or len(active) == len(games):
            return margin, games, edges, stabilizers, rows
        candidate = restrict_family(games, edges, stabilizers, active)
        candidate_margin = solve_quotient(*candidate, 6)
        if abs(candidate_margin - margin) > 1e-8 * max(1.0, abs(margin)):
            rows[-1]["rejected_candidate_margin"] = candidate_margin
            return margin, games, edges, stabilizers, rows
        games, edges, stabilizers = candidate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--mode", choices=("null", "veto"), default="veto")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source_games, source_edges, _source_stabilizers = load_quotient(args.input)
    lift = veto_lift if args.mode == "veto" else dummy_lift
    representatives = [lift(game) for game in source_games]
    permutations = list(itertools.permutations(range(6)))
    orbit_games = []
    orbit_maps = []
    orbit_index = {}
    stabilizers = [[] for _ in representatives]
    for game_index, game in enumerate(representatives):
        for permutation in permutations:
            image = permute_game(game, permutation)
            if image == game:
                stabilizers[game_index].append(permutation)
            if image in orbit_index:
                continue
            orbit_index[image] = len(orbit_games)
            orbit_games.append(image)
            inverse = [0] * 6
            for source, target in enumerate(permutation):
                inverse[target] = source
            orbit_maps.append((game_index, tuple(inverse)))
    print(
        json.dumps(
            {
                "phase": "orbit",
                "representatives": len(representatives),
                "orbit_games": len(orbit_games),
            }
        ),
        flush=True,
    )
    raw_edges = add_induced_edges(orbit_games)
    quotient_edges = set()
    for lower, upper, coalition in raw_edges:
        lower_game, lower_inverse = orbit_maps[lower]
        upper_game, upper_inverse = orbit_maps[upper]
        for player in range(6):
            if coalition >> player & 1:
                quotient_edges.add(
                    (
                        lower_game,
                        lower_inverse[player],
                        upper_game,
                        upper_inverse[player],
                    )
                )
    print(
        json.dumps(
            {
                "phase": "comparisons",
                "raw_edges": len(raw_edges),
                "quotient_edges": len(quotient_edges),
            }
        ),
        flush=True,
    )

    if not quotient_edges:
        output = {
            "status": "no_comparisons_survive_dummy_lift",
            "source": str(args.input),
            "mode": args.mode,
            "source_games": len(source_games),
            "source_quotient_edges": len(source_edges),
            "orbit_games": len(orbit_games),
            "raw_edges": len(raw_edges),
            "quotient_edges": 0,
        }
        args.output.write_text(json.dumps(output, indent=2) + "\n")
        print(json.dumps(output, indent=2))
        return 0

    best = None
    for nodes in components(len(representatives), quotient_edges):
        local = restrict_family(
            representatives, quotient_edges, stabilizers, nodes
        )
        margin = solve_quotient(*local, 6)
        if best is None or margin < best[0]:
            best = (margin, *local)
        if margin < -1e-8:
            break
    if best is None:
        raise RuntimeError("six-player lift has no component")
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
        "mode": args.mode,
        "source_games": len(source_games),
        "orbit_games": len(orbit_games),
        "raw_edges": len(raw_edges),
        "quotient_edges": len(quotient_edges),
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
