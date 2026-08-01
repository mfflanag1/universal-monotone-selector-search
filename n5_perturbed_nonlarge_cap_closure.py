#!/usr/bin/env python3
"""Grow an exact boundary closure above a perturbed non-large aggregate cap."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from n5_facet_search import load_facets
from n5_sparse_family_lp import sparse_family_slack_float
from nonlarge_complex_search import biswas_nonlarge_game


F = Fraction
N = 5
GRAND = 31


def game_key(game: np.ndarray) -> tuple[float, ...]:
    return tuple(float(round(value, 10)) for value in game)


def moves(game: np.ndarray, facets: np.ndarray):
    slack = facets @ game
    for coalition in range(1, GRAND):
        negative = facets[:, coalition] < 0
        bounds = [
            float(np.min(slack[negative] / -facets[negative, coalition]))
        ] if negative.any() else []
        bounds.extend(
            game[upper] - game[coalition]
            for upper in range(1, 32)
            if upper != coalition and upper & coalition == coalition
        )
        delta = min(bounds)
        if delta > 1e-9:
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
        if delta > 1e-9:
            neighbor = game.copy()
            neighbor[coalition] -= delta
            yield neighbor, coalition, False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epsilon", type=F, default=F(1, 100))
    parser.add_argument("--grand-bump", type=F, default=F(1, 10000))
    parser.add_argument("--depth", type=int, default=2)
    parser.add_argument("--max-nodes", type=int, default=5000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facets = np.asarray(load_facets()[0], dtype=float)
    base = biswas_nonlarge_game(F(2))
    convex = [F(2 ** mask.bit_count() - mask.bit_count() - 1) for mask in range(32)]
    lower = np.asarray(
        [float(base[mask] + args.epsilon * convex[mask]) for mask in range(32)]
    )
    upper = lower.copy()
    upper[GRAND] += float(args.grand_bump)
    games = [lower, upper]
    indices = {game_key(lower): 0, game_key(upper): 1}
    edges: set[tuple[int, int, int]] = {(0, 1, GRAND)}
    frontier = [1]
    depth_rows = []
    for depth in range(1, args.depth + 1):
        next_frontier = []
        for node in frontier:
            for neighbor, coalition, increase in moves(games[node], facets):
                neighbor_key = game_key(neighbor)
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
        depth_rows.append(
            {"depth": depth, "nodes": len(games), "edges": len(edges)}
        )
        print(json.dumps(depth_rows[-1]), flush=True)
        if not frontier or len(games) >= args.max_nodes:
            break
    margin = sparse_family_slack_float(
        games, sorted(edges), N, "highs-ipm"
    )
    if margin is None:
        raise RuntimeError("sparse family LP failed")
    payload = {
        "status": "perturbed_nonlarge_cap_closure_float",
        "epsilon": str(args.epsilon),
        "grand_bump": str(args.grand_bump),
        "depth_rows": depth_rows,
        "nodes": len(games),
        "edges": len(edges),
        "margin": margin,
        "games": [[float(value) for value in game] for game in games],
        "edge_list": [list(edge) for edge in sorted(edges)],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"margin": margin, "nodes": len(games), "edges": len(edges)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
