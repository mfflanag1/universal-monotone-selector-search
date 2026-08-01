#!/usr/bin/env python3
"""Explore legal point-deletion networks from a full permutation orbit."""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np

sys.path.append(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)

from n5_point_deletion_family_search import connected_components
from n5_sparse_family_lp import sparse_family_slack_float


def coalition_values(point: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(value for player, value in enumerate(point) if coalition >> player & 1)
        for coalition in range(1 << len(point))
    )


def minima_and_unique(point_values, active, coalition_count):
    active_values = point_values[active]
    minima = active_values.min(axis=0)
    minimizing = active_values == minima
    counts = minimizing.sum(axis=0)
    unique = [[] for _ in point_values]
    uniquely_minimized = np.flatnonzero(counts == 1)
    for coalition in uniquely_minimized:
        if 0 < coalition < coalition_count - 1:
            local_point = int(np.flatnonzero(minimizing[:, coalition])[0])
            unique[active[local_point]].append(int(coalition))
    return tuple(int(value) for value in minima), unique


def add_induced_edges(games, edges, coalition_count):
    for coalition in range(1, coalition_count - 1):
        groups = {}
        for index, game in enumerate(games):
            signature = game[:coalition] + game[coalition + 1 :]
            groups.setdefault(signature, []).append(index)
        for indices in groups.values():
            by_value = {}
            for index in indices:
                by_value.setdefault(games[index][coalition], index)
            ordered = [by_value[value] for value in sorted(by_value)]
            for lower_position, lower in enumerate(ordered):
                for upper in ordered[lower_position + 1 :]:
                    edges.add((lower, upper, coalition))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", required=True)
    parser.add_argument("--walks", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--solve-every", type=int, default=100)
    parser.add_argument("--restart-probability", type=float, default=0.35)
    parser.add_argument("--state-reservoir", type=int, default=100_000)
    parser.add_argument("--cache-size", type=int, default=50_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pattern = tuple(int(value) for value in args.pattern.split(","))
    n = len(pattern)
    coalition_count = 1 << n
    points = sorted(
        {
            tuple(pattern[index] for index in permutation)
            for permutation in itertools.permutations(range(n))
        }
    )
    point_values = np.asarray(
        [coalition_values(point) for point in points], dtype=np.int16
    )
    full = (1 << len(points)) - 1
    rng = random.Random(args.seed)
    accepted_states = [full]
    state_evaluations = 0
    games = []
    game_index = {}
    edges = set()
    best = None
    state_cache = OrderedDict()

    def state_info(state):
        nonlocal state_evaluations
        cached = state_cache.get(state)
        if cached is not None:
            state_cache.move_to_end(state)
            return cached
        state_evaluations += 1
        active = [index for index in range(len(points)) if state >> index & 1]
        game, unique = minima_and_unique(point_values, active, coalition_count)
        legal = [point for point in active if len(unique[point]) <= 1]
        productive = [point for point in legal if unique[point]]
        coalitions = {
            point: unique[point][0]
            for point in productive
        }
        cached = (game, legal, productive, coalitions)
        state_cache[state] = cached
        if len(state_cache) > args.cache_size:
            state_cache.popitem(last=False)
        return cached

    for walk in range(1, args.walks + 1):
        state = full if rng.random() < args.restart_probability else rng.choice(accepted_states)
        while state.bit_count() > 1:
            before, legal, productive, coalitions = state_info(state)
            if not legal:
                break
            removed = rng.choice(productive if productive and rng.random() < 0.8 else legal)
            successor = state ^ (1 << removed)
            after = state_info(successor)[0]
            if len(accepted_states) < args.state_reservoir:
                accepted_states.append(successor)
            else:
                replacement = rng.randrange(state_evaluations + 1)
                if replacement < args.state_reservoir:
                    accepted_states[replacement] = successor
            if removed in coalitions:
                coalition = coalitions[removed]
                for game in (before, after):
                    if game not in game_index:
                        game_index[game] = len(games)
                        games.append(game)
                edges.add((game_index[before], game_index[after], coalition))
            state = successor

        if walk % args.solve_every and walk != args.walks:
            continue
        add_induced_edges(games, edges, coalition_count)
        components = connected_components(len(games), edges)
        current = None
        for nodes, local_edges in components[:20]:
            margin = sparse_family_slack_float(
                [games[node] for node in nodes], local_edges, n, "highs-ipm"
            )
            if margin is not None and (current is None or margin < current):
                current = margin
            if margin is not None and (best is None or margin < best[0]):
                best = (
                    margin,
                    [games[node] for node in nodes],
                    local_edges,
                )
            if margin is not None and margin < -1e-8:
                break
        print(
            json.dumps(
                {
                    "walk": walk,
                    "state_evaluations": state_evaluations,
                    "state_reservoir": len(accepted_states),
                    "games": len(games),
                    "edges": len(edges),
                    "components": len(components),
                    "largest_component": len(components[0][0]) if components else 0,
                    "margin": current,
                }
            ),
            flush=True,
        )
        if current is not None and current < -1e-8:
            break

    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "pattern": list(pattern),
        "orbit_size": len(points),
        "state_evaluations": state_evaluations,
        "state_reservoir_size": len(accepted_states),
        "game_count": len(games),
        "edge_count": len(edges),
        "best": None
        if best is None
        else {
            "margin": best[0],
            "games": [list(game) for game in best[1]],
            "edges": [list(edge) for edge in best[2]],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "orbit_size", "state_evaluations", "state_reservoir_size", "game_count", "edge_count")}, indent=2))
    return 1 if best is not None and best[0] < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
