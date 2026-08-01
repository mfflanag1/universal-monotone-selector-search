#!/usr/bin/env python3
"""Explore the full 30-point BPPV orbit by legal or invisible deletions."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from fractions import Fraction
from pathlib import Path

from n5_point_deletion_family_search import connected_components
from n5_sparse_family_lp import sparse_family_slack_float


F = Fraction
N = 5
GRAND = (1 << N) - 1


def orbit_points():
    point = (3, 1, 1, 2, 2)
    return sorted(
        {
            tuple(F(point[index]) for index in permutation)
            for permutation in itertools.permutations(range(N))
        }
    )


def state_game(point_values, state, cache):
    if state not in cache:
        active = [
            index
            for index in range(len(point_values))
            if state >> index & 1
        ]
        cache[state] = tuple(
            min(point_values[index][coalition] for index in active)
            for coalition in range(GRAND + 1)
        )
    return cache[state]


def add_induced_edges(games, edges):
    for coalition in range(1, GRAND):
        groups = {}
        for index, game in enumerate(games):
            signature = game[:coalition] + game[coalition + 1 :]
            groups.setdefault(signature, []).append(index)
        for indices in groups.values():
            ordered = sorted(indices, key=lambda index: games[index][coalition])
            for lower_position, lower in enumerate(ordered):
                for upper in ordered[lower_position + 1 :]:
                    edges.add((lower, upper, coalition))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--walks", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--solve-every", type=int, default=50)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    points = orbit_points()
    point_values = [
        tuple(
            sum(
                point[player]
                for player in range(N)
                if coalition >> player & 1
            )
            for coalition in range(GRAND + 1)
        )
        for point in points
    ]
    full = (1 << len(points)) - 1
    cache = {}
    games = []
    game_index = {}
    edges = set()
    accepted_states = [full]
    best_margin = None
    for walk in range(1, args.walks + 1):
        state = full if rng.random() < 0.7 else rng.choice(accepted_states)
        while state.bit_count() > 1:
            before = state_game(point_values, state, cache)
            candidates = [
                index for index in range(len(points)) if state >> index & 1
            ]
            rng.shuffle(candidates)
            moves = []
            for removed in candidates:
                successor = state ^ (1 << removed)
                after = state_game(point_values, successor, cache)
                differences = [
                    coalition
                    for coalition, (left, right) in enumerate(
                        zip(before, after, strict=True)
                    )
                    if left != right
                ]
                if not differences:
                    moves.append((0, successor, after))
                elif (
                    len(differences) == 1
                    and differences[0] not in (0, GRAND)
                    and before[differences[0]] < after[differences[0]]
                ):
                    moves.append((differences[0], successor, after))
            if not moves:
                break
            coalition, successor, after = rng.choice(moves)
            accepted_states.append(successor)
            if coalition:
                for game in (before, after):
                    if game not in game_index:
                        game_index[game] = len(games)
                        games.append(game)
                edges.add(
                    (game_index[before], game_index[after], coalition)
                )
            state = successor
        if walk % args.solve_every and walk != args.walks:
            continue
        add_induced_edges(games, edges)
        components = connected_components(len(games), edges)
        best_margin = None
        for nodes, local_edges in components:
            margin = sparse_family_slack_float(
                [games[node] for node in nodes], local_edges, N
            )
            if margin is not None and (
                best_margin is None or margin < best_margin
            ):
                best_margin = margin
            if margin is not None and margin < -1e-8:
                break
        print(
            json.dumps(
                {
                    "walk": walk,
                    "states": len(cache),
                    "games": len(games),
                    "edges": len(edges),
                    "components": len(components),
                    "margin": best_margin,
                }
            ),
            flush=True,
        )
        if best_margin is not None and best_margin < -1e-8:
            break
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best_margin is not None and best_margin < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "point_count": len(points),
        "explored_state_count": len(cache),
        "game_count": len(games),
        "edge_count": len(edges),
        "best_margin": best_margin,
        "games": [[str(value) for value in game] for game in games],
        "edges": [list(edge) for edge in sorted(edges)],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
