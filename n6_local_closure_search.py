#!/usr/bin/env python3
"""Grow depth-two exact-coordinate closures around six-player point games."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np

from family_search import family_slack_float
from n6_local_star_search import GRAND, compositions, point_game
from n6_mixed_terminal_search import load_facets


def key(game: np.ndarray) -> tuple[float, ...]:
    return tuple(float(round(value, 9)) for value in game)


def moves(game: np.ndarray, facets: np.ndarray):
    slack = facets @ game
    for coalition in range(1, GRAND):
        negative = facets[:, coalition] < 0
        bounds = [
            float(np.min(slack[negative] / -facets[negative, coalition]))
        ] if negative.any() else []
        bounds.extend(
            game[upper] - game[coalition]
            for upper in range(1, 64)
            if upper != coalition and upper & coalition == coalition
        )
        delta = min(bounds)
        if delta > 1e-8:
            neighbor = game.copy()
            neighbor[coalition] += delta
            yield neighbor, coalition, True

        positive = facets[:, coalition] > 0
        bounds = [
            float(np.min(slack[positive] / facets[positive, coalition]))
        ] if positive.any() else []
        bounds.extend(
            game[coalition] - game[lower]
            for lower in range(GRAND)
            if lower != coalition and lower & coalition == lower
        )
        delta = min(bounds)
        if delta > 1e-8:
            neighbor = game.copy()
            neighbor[coalition] -= delta
            yield neighbor, coalition, False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=10)
    parser.add_argument("--depth", type=int, default=2)
    parser.add_argument("--point-count", type=int, default=10)
    parser.add_argument("--total", type=int, default=10)
    parser.add_argument("--max-nodes", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pool = list(compositions(args.total))
    facets = np.asarray(load_facets(), dtype=float)
    rows = []
    best = None
    for trial in range(args.trials):
        root = point_game(rng.sample(pool, args.point_count))
        games = [root]
        indices = {key(root): 0}
        edges: set[tuple[int, int, int]] = set()
        frontier = [0]
        for _depth in range(args.depth):
            next_frontier = []
            for node in frontier:
                for neighbor, coalition, increase in moves(games[node], facets):
                    neighbor_key = key(neighbor)
                    target = indices.get(neighbor_key)
                    if target is None:
                        if len(games) >= args.max_nodes:
                            continue
                        target = len(games)
                        indices[neighbor_key] = target
                        games.append(neighbor)
                        next_frontier.append(target)
                    edge = (
                        (node, target, coalition)
                        if increase
                        else (target, node, coalition)
                    )
                    if edge[0] != edge[1]:
                        edges.add(edge)
            frontier = next_frontier
            if not frontier or len(games) >= args.max_nodes:
                break
        margin, _ = family_slack_float(games, sorted(edges), 6)
        row = {
            "trial": trial,
            "nodes": len(games),
            "edges": len(edges),
            "margin": margin,
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
        if best is None or margin < best[0]:
            best = (margin, games, sorted(edges))
        if margin < -1e-8:
            break
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "rows": rows,
        "best": None
        if best is None
        else {
            "margin": best[0],
            "games": [[float(value) for value in game] for game in best[1]],
            "edges": [list(edge) for edge in best[2]],
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
