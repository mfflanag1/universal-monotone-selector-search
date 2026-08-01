#!/usr/bin/env python3
"""Search legal coordinate paths between cloned exact Housman--Clark cores."""

from __future__ import annotations

import argparse
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_sparse_family_lp import sparse_family_slack_float


def core_support(game: list[float], objective: np.ndarray) -> float:
    n = 4
    grand = 15
    rows = []
    rhs = []
    for coalition in range(1, grand):
        rows.append(
            [-float(bool(coalition >> player & 1)) for player in range(n)]
        )
        rhs.append(-game[coalition])
    result = linprog(
        objective,
        A_ub=np.asarray(rows),
        b_ub=np.asarray(rhs),
        A_eq=np.ones((1, n)),
        b_eq=[game[grand]],
        bounds=[(None, None)] * n,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return float(result.fun)


def cloned_exact_game(game: list[float], pair_totals: tuple[float, ...]) -> list[float]:
    n = 8
    result = []
    for coalition in range(1 << n):
        objective = np.zeros(4)
        constant = 0.0
        for player in range(4):
            if coalition >> player & 1:
                objective[player] += 1.0
            if coalition >> (player + 4) & 1:
                objective[player] -= 1.0
                constant += pair_totals[player]
        result.append(constant + core_support(game, objective))
    return result


def interior_game(n: int) -> list[float]:
    return [coalition.bit_count() ** 2 / n**2 for coalition in range(1 << n)]


def orderings(
    changed: list[int], active: int, rng: random.Random, variants: int
) -> list[list[int]]:
    originals = active
    clones = active << 4

    def score(mask: int) -> tuple[int, int, int]:
        good = (mask & originals).bit_count()
        bad = (mask & clones).bit_count()
        return good - bad, good, -bad

    candidates = [
        sorted(changed, key=score),
        sorted(changed, key=score, reverse=True),
    ]
    while len(candidates) < variants:
        candidate = list(changed)
        rng.shuffle(candidate)
        candidates.append(candidate)
    return candidates[:variants]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--rounds", type=int, default=1)
    parser.add_argument("--epsilon", type=float, default=0.0)
    parser.add_argument("--buffer", choices=("convex", "grand"), default="convex")
    parser.add_argument("--variants", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--save-family", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = json.loads(args.source.read_text())["family"]
    old_games = [
        [float(Fraction(value)) for value in game] for game in raw["games"]
    ]
    pair_totals = (8.0, 8.0, 8.0, 8.0)
    endpoints = [cloned_exact_game(game, pair_totals) for game in old_games]
    h = interior_game(8) if args.buffer == "convex" else [0.0] * 255 + [1.0]
    endpoints = [
        [value + args.epsilon * h[coalition] for coalition, value in enumerate(game)]
        for game in endpoints
    ]
    games = [endpoints[0]]
    edges: list[tuple[int, int, int]] = []
    rng = random.Random(args.seed)
    terminal_nodes = []
    for branch, edge in enumerate(raw["edges"], start=1):
        upper_game = endpoints[int(edge["upper"])]
        delta = [
            upper_game[coalition] - endpoints[0][coalition]
            for coalition in range(256)
        ]
        changed = [
            coalition
            for coalition in range(1, 255)
            if delta[coalition] > 1e-9
        ]
        active = int(edge["coalition"])
        branch_terminals = []
        for ordering in orderings(changed, active, rng, args.variants):
            current_game = list(endpoints[0])
            current_node = 0
            for _round in range(args.rounds):
                for coalition in ordering:
                    next_game = list(current_game)
                    next_game[coalition] += delta[coalition] / args.rounds
                    if (
                        _round == args.rounds - 1
                        and coalition == ordering[-1]
                        and not branch_terminals
                    ):
                        next_node = len(games)
                        games.append(next_game)
                        branch_terminals.append(next_node)
                    elif (
                        _round == args.rounds - 1
                        and coalition == ordering[-1]
                        and branch_terminals
                    ):
                        next_node = branch_terminals[0]
                    else:
                        next_node = len(games)
                        games.append(next_game)
                    edges.append((current_node, next_node, coalition))
                    current_game = next_game
                    current_node = next_node
        terminal_nodes.append(branch_terminals[0])
    margin = sparse_family_slack_float(games, edges, 8, "highs-ipm")
    payload = {
        "status": (
            "incompatible_cloned_path_family_found"
            if margin is not None and margin < -1e-8
            else "cloned_path_family_compatible"
        ),
        "rounds": args.rounds,
        "epsilon": args.epsilon,
        "buffer": args.buffer,
        "variants": args.variants,
        "game_count": len(games),
        "edge_count": len(edges),
        "changed_coordinates": [
            sum(
                endpoints[branch][coalition] - endpoints[0][coalition] > 1e-9
                for coalition in range(1, 255)
            )
            for branch in range(1, len(endpoints))
        ],
        "margin": margin,
        "terminal_nodes": terminal_nodes,
    }
    if args.save_family:
        payload["games"] = games
        payload["edges"] = [list(edge) for edge in edges]
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {key: value for key, value in payload.items() if key not in {"games", "edges"}},
            indent=2,
        )
    )
    return 1 if margin is not None and margin < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
