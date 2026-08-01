#!/usr/bin/env python3
"""Explore the exact {0,1,2,3}-grid components reachable from grand-2 games."""

from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path

import numpy as np

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_boolean_exact_domain import GRAND, N
from n5_facet_search import load_facets
from n5_point_deletion_family_search import connected_components
from n5_sparse_family_lp import sparse_family_slack_float
from n5_symmetric_family_lp import symmetric_family_result
from n5_ternary_exact_domain import monotone_boolean_codes, vector


def exact_grand_two_games(facet_matrix: np.ndarray) -> set[bytes]:
    codes = monotone_boolean_codes()
    boolean_games = np.stack([vector(code) for code in codes])
    scores = boolean_games[:, 1:].astype(np.int16) @ facet_matrix.T
    games: set[bytes] = {bytes(GRAND + 1)}
    for index in np.flatnonzero(np.all(scores >= 0, axis=1)):
        games.add(bytes(boolean_games[index]))
    for upper_index, upper_code in enumerate(codes):
        nested = np.fromiter(
            (index for index, lower_code in enumerate(codes) if lower_code & ~upper_code == 0),
            dtype=np.int32,
        )
        feasible = nested[np.all(scores[nested] + scores[upper_index] >= 0, axis=1)]
        upper = boolean_games[upper_index]
        for lower_index in feasible:
            games.add(bytes(upper + boolean_games[lower_index]))
    return games


def local_monotone(game: bytes, coalition: int, value: int) -> bool:
    for player in range(N):
        bit = 1 << player
        if coalition & bit:
            if game[coalition ^ bit] > value:
                return False
        elif game[coalition | bit] < value:
            return False
    return True


def reachable_section(
    facet_matrix: np.ndarray,
    seeds: set[bytes],
    grand_worth: int,
    max_games: int,
) -> tuple[list[bytes], bool]:
    games: list[bytes] = []
    scores: list[np.ndarray] = []
    index: dict[bytes, int] = {}
    queue: deque[int] = deque()
    for seed in sorted(seeds):
        raised = bytearray(seed)
        raised[GRAND] = grand_worth
        game = bytes(raised)
        if game in index:
            continue
        score = facet_matrix @ np.frombuffer(game, dtype=np.uint8)[1:].astype(np.int16)
        if np.any(score < 0):
            raise AssertionError("grand-worth raise unexpectedly left the exact cone")
        index[game] = len(games)
        games.append(game)
        scores.append(score)
        queue.append(len(games) - 1)

    complete = True
    expanded = 0
    while queue:
        game_index = queue.popleft()
        game = games[game_index]
        score = scores[game_index]
        for coalition in range(1, GRAND):
            current = game[coalition]
            for delta in (-1, 1):
                value = current + delta
                if value < 0 or value > grand_worth or not local_monotone(game, coalition, value):
                    continue
                candidate_score = score + delta * facet_matrix[:, coalition - 1]
                if np.any(candidate_score < 0):
                    continue
                candidate_array = bytearray(game)
                candidate_array[coalition] = value
                candidate = bytes(candidate_array)
                if candidate in index:
                    continue
                if len(games) >= max_games:
                    complete = False
                    queue.clear()
                    break
                index[candidate] = len(games)
                games.append(candidate)
                scores.append(candidate_score)
                queue.append(len(games) - 1)
            if not complete:
                break
        expanded += 1
        if expanded % 10000 == 0:
            print(
                json.dumps(
                    {"phase": "bfs", "expanded": expanded, "games": len(games), "queue": len(queue)}
                ),
                flush=True,
            )
    return games, complete


def induced_edges(
    games: list[bytes], max_value: int, include_grand: bool = False
) -> set[tuple[int, int, int]]:
    index = {game: position for position, game in enumerate(games)}
    edges: set[tuple[int, int, int]] = set()
    for lower, game in enumerate(games):
        for coalition in range(1, GRAND + int(include_grand)):
            for value in range(game[coalition] + 1, max_value + 1):
                successor = bytearray(game)
                successor[coalition] = value
                upper = index.get(bytes(successor))
                if upper is not None:
                    edges.add((lower, upper, coalition))
    return edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-games", type=int, default=500_000)
    parser.add_argument("--solve-max-games", type=int, default=150_000)
    parser.add_argument("--cap", type=int, default=3)
    parser.add_argument("--symmetric-solve", action="store_true")
    parser.add_argument("--include-lower-sections", action="store_true")
    args = parser.parse_args()

    facets, _types = load_facets()
    facet_matrix = np.asarray(facets, dtype=np.int16)[:, 1:]
    seeds = exact_grand_two_games(facet_matrix)
    if args.cap < 3 or args.cap > 255:
        raise ValueError("cap must be between 3 and 255")
    lower_sections = sorted(seeds)
    section_seeds = seeds
    games = []
    complete = True
    for grand_worth in range(3, args.cap + 1):
        games, complete = reachable_section(
            facet_matrix, section_seeds, grand_worth, args.max_games
        )
        print(
            json.dumps(
                {
                    "phase": "section_complete" if complete else "section_truncated",
                    "grand_worth": grand_worth,
                    "games": len(games),
                }
            ),
            flush=True,
        )
        if not complete or grand_worth == args.cap:
            break
        lower_sections.extend(games)
        section_seeds = set(games)
    if not complete:
        payload = {
            "status": "reachable_section_truncated",
            "scope": f"grand-{grand_worth} reachable exact integer grid",
            "grand_two_seed_count": len(seeds),
            "cap": args.cap,
            "max_games": args.max_games,
            "game_count_at_truncation": len(games),
        }
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
        print(json.dumps(payload, indent=2))
        return 2
    if args.include_lower_sections:
        games = lower_sections + games
    edges = induced_edges(games, args.cap, args.include_lower_sections)
    if args.symmetric_solve:
        if not complete:
            raise RuntimeError("symmetry quotient requires a complete permutation-closed BFS")
        result = symmetric_family_result(games, edges)
        if not result.success:
            raise RuntimeError(result.message)
        margin = float(result.x[-1])
        payload = {
            "status": (
                "incompatible_exact_family_found"
                if margin < -1e-8
                else "no_incompatible_exact_family_found"
            ),
            "scope": (
                f"combined grand-at-most-{args.cap} reachable exact integer grid"
                if args.include_lower_sections
                else f"grand-{args.cap} exact integer-grid components reachable from lower sections"
            ),
            "grand_two_seed_count": len(seeds),
            "bfs_complete": complete,
            "game_count": len(games),
            "legal_edge_count": len(edges),
            "symmetry_quotient": result.quotient_counts,
            "margin": margin,
        }
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
        print(json.dumps(payload, indent=2))
        return 1 if margin < -1e-8 else 0
    components = connected_components(len(games), edges)
    rows = []
    best = None
    for component, (nodes, local_edges) in enumerate(components):
        margin = None
        if len(nodes) <= args.solve_max_games:
            local_games = [tuple(games[node]) for node in nodes]
            margin = sparse_family_slack_float(local_games, local_edges, N, "highs-ipm")
            if margin is not None and (best is None or margin < best[0]):
                best = (margin, local_games, local_edges)
        row = {
            "component": component,
            "nodes": len(nodes),
            "edges": len(local_edges),
            "cycle_rank": len(local_edges) - len(nodes) + 1,
            "margin": margin,
        }
        rows.append(row)
        print(json.dumps(row), flush=True)

    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "scope": "grand-3 exact integer-grid components reachable from all grand-at-most-2 games",
        "grand_two_seed_count": len(seeds),
        "bfs_complete": complete,
        "game_count": len(games),
        "legal_edge_count": len(edges),
        "components": rows,
        "best": None
        if best is None
        else {
            "margin": best[0],
            "games": [list(game) for game in best[1]] if best[0] < -1e-8 else None,
            "edges": [list(edge) for edge in best[2]] if best[0] < -1e-8 else None,
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if best is not None and best[0] < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
