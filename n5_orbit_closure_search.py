#!/usr/bin/env python3
"""Search induced comparison closures of player-permuted n=5 gadgets."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

from n5_facet_search import exact_archive, permute_mask, search


F = Fraction
N = 5
GRAND = (1 << N) - 1


def parse_players(text: str) -> tuple[int, ...]:
    players = tuple(int(value) for value in text.split(",") if value)
    if not players or len(set(players)) != len(players):
        raise ValueError("permuted players must be distinct and nonempty")
    if any(player not in range(N) for player in players):
        raise ValueError("players must be in 0,1,2,3,4")
    return players


def subgroup(players: tuple[int, ...]) -> list[tuple[int, ...]]:
    permutations = []
    for image in itertools.permutations(players):
        permutation = list(range(N))
        for old, new in zip(players, image, strict=True):
            permutation[old] = new
        permutations.append(tuple(permutation))
    return permutations


def normalized_games(payload: dict) -> list[tuple[Fraction, ...]]:
    raw = payload["family"]
    games = [
        tuple(F(value) for value in game)
        for game in raw["games"]
    ]
    normalizer = games[0][GRAND]
    if normalizer <= 0:
        raise ValueError("grand-worth normalizer must be positive")
    return [
        tuple(value / normalizer for value in game)
        for game in games
    ]


def permuted_displacement(
    game: tuple[Fraction, ...],
    anchor: tuple[Fraction, ...],
    permutation: tuple[int, ...],
) -> tuple[Fraction, ...]:
    displacement = [F(0)] * (1 << N)
    for coalition in range(1 << N):
        displacement[permute_mask(coalition, permutation)] = (
            game[coalition] - anchor[coalition]
        )
    return tuple(displacement)


def induced_topology(
    games: list[tuple[Fraction, ...]],
    anchor: int,
    permutations: list[tuple[int, ...]],
) -> tuple[list[tuple[int, int, int]], list[float], int]:
    anchor_game = games[anchor]
    state_ids: dict[tuple[Fraction, ...], int] = {}
    for permutation in permutations:
        for game in games:
            state = permuted_displacement(
                game, anchor_game, permutation
            )
            state_ids.setdefault(state, len(state_ids))

    states: list[tuple[Fraction, ...] | None] = [None] * len(state_ids)
    for state, index in state_ids.items():
        states[index] = state
    exact_states = [
        state for state in states if state is not None
    ]

    edges: list[tuple[int, int, int]] = []
    bumps: list[float] = []
    for left in range(len(exact_states)):
        for right in range(left + 1, len(exact_states)):
            changed = [
                coalition
                for coalition in range(1, GRAND + 1)
                if exact_states[left][coalition]
                != exact_states[right][coalition]
            ]
            if len(changed) != 1:
                continue
            coalition = changed[0]
            delta = (
                exact_states[right][coalition]
                - exact_states[left][coalition]
            )
            if delta > 0:
                edges.append((left, right, coalition))
                bumps.append(float(delta))
            else:
                edges.append((right, left, coalition))
                bumps.append(float(-delta))
    if not edges:
        raise RuntimeError("orbit closure has no comparison edges")
    zero = tuple(F(0) for _ in range(1 << N))
    root = state_ids[zero]
    adjacency: dict[int, set[int]] = {
        node: set() for node in range(len(exact_states))
    }
    for lower, upper, _ in edges:
        adjacency[lower].add(upper)
        adjacency[upper].add(lower)
    component = {root}
    frontier = [root]
    while frontier:
        node = frontier.pop()
        for neighbor in adjacency[node] - component:
            component.add(neighbor)
            frontier.append(neighbor)
    remap = {
        old: new for new, old in enumerate(sorted(component))
    }
    connected_edges = []
    connected_bumps = []
    for edge, bump in zip(edges, bumps, strict=True):
        lower, upper, coalition = edge
        if lower in component and upper in component:
            connected_edges.append(
                (remap[lower], remap[upper], coalition)
            )
            connected_bumps.append(bump)
    return connected_edges, connected_bumps, len(component)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--anchor", type=int, required=True)
    parser.add_argument("--permuted-players", required=True)
    parser.add_argument("--starts", type=int, default=6)
    parser.add_argument("--iterations", type=int, default=18)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--max-denominator", type=int, default=1_000_000)
    parser.add_argument("--float-only", action="store_true")
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    games = normalized_games(payload)
    if args.anchor not in range(len(games)):
        raise ValueError("anchor is outside the gadget")
    players = parse_players(args.permuted_players)
    permutations = subgroup(players)
    edges, bumps, node_count = induced_topology(
        games, args.anchor, permutations
    )
    print(
        json.dumps(
            {
                "node_count": node_count,
                "edge_count": len(edges),
                "permutation_count": len(permutations),
            }
        ),
        flush=True,
    )
    record = search(
        edges,
        bumps,
        args.starts,
        args.iterations,
        args.seed,
    )
    record["orbit_closure"] = {
        "source": str(args.input),
        "anchor": args.anchor,
        "permuted_players": list(players),
        "permutations": [
            list(permutation) for permutation in permutations
        ],
    }
    if args.float_only:
        result = {
            "status": "float_orbit_closure",
            "record": record,
        }
    else:
        result = exact_archive(
            record,
            "n5_induced_player_orbit_closure",
            args.max_denominator,
        )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "node_count": node_count,
                "edge_count": len(edges),
                "best_margin": (
                    record["best_margin_float"]
                    if args.float_only
                    else result[
                        "max_min_monotonicity_margin_exact"
                    ]
                ),
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
