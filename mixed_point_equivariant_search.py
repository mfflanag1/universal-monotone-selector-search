#!/usr/bin/env python3
"""Search permutation-equivariant obstructions in exact point-deletion families."""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from mixed_point_family_fast import build, random_points
from n5_intrinsic_permutation_union import permute_game, permute_mask
from n5_point_deletion_family_search import connected_components


def equivariant_margin(
    base_games: list[tuple[int, ...]],
    base_edges: list[tuple[int, int, int]],
    n: int,
) -> tuple[float, int, int]:
    grand = (1 << n) - 1
    permutations = list(itertools.permutations(range(n)))
    games: list[tuple[int, ...]] = []
    game_index: dict[tuple[int, ...], int] = {}
    images: dict[tuple[int, tuple[int, ...]], int] = {}
    for base_index, game in enumerate(base_games):
        for permutation in permutations:
            image = tuple(int(value) for value in permute_game(game, permutation))
            if image not in game_index:
                game_index[image] = len(games)
                games.append(image)
            images[base_index, permutation] = game_index[image]
    edges = {
        (
            images[lower, permutation],
            images[upper, permutation],
            permute_mask(coalition, permutation),
        )
        for permutation in permutations
        for lower, upper, coalition in base_edges
    }
    allocation_count = len(games) * n
    slack = allocation_count
    ub_rows: list[int] = []
    ub_columns: list[int] = []
    ub_values: list[float] = []
    ub_rhs: list[float] = []
    row = 0
    for game_index_value, game in enumerate(games):
        for coalition in range(1, grand):
            for player in range(n):
                if coalition >> player & 1:
                    ub_rows.append(row)
                    ub_columns.append(game_index_value * n + player)
                    ub_values.append(-1.0)
            ub_rhs.append(-float(game[coalition]))
            row += 1
    for lower, upper, coalition in sorted(edges):
        for player in range(n):
            if not coalition >> player & 1:
                continue
            ub_rows.extend((row, row, row))
            ub_columns.extend((lower * n + player, upper * n + player, slack))
            ub_values.extend((1.0, -1.0, 1.0))
            ub_rhs.append(0.0)
            row += 1
    inequalities = coo_matrix(
        (ub_values, (ub_rows, ub_columns)),
        shape=(row, allocation_count + 1),
    ).tocsr()

    eq_rows: list[int] = []
    eq_columns: list[int] = []
    eq_values: list[float] = []
    eq_rhs: list[float] = []
    row = 0
    for game_index_value, game in enumerate(games):
        for player in range(n):
            eq_rows.append(row)
            eq_columns.append(game_index_value * n + player)
            eq_values.append(1.0)
        eq_rhs.append(float(game[grand]))
        row += 1
    identity = tuple(range(n))
    for base_index in range(len(base_games)):
        base_image = images[base_index, identity]
        for permutation in permutations:
            image = images[base_index, permutation]
            for player in range(n):
                eq_rows.extend((row, row))
                eq_columns.extend(
                    (image * n + permutation[player], base_image * n + player)
                )
                eq_values.extend((1.0, -1.0))
                eq_rhs.append(0.0)
                row += 1
    equalities = coo_matrix(
        (eq_values, (eq_rows, eq_columns)),
        shape=(row, allocation_count + 1),
    ).tocsr()
    objective = np.zeros(allocation_count + 1)
    objective[slack] = -1.0
    result = linprog(
        objective,
        A_ub=inequalities,
        b_ub=np.asarray(ub_rhs),
        A_eq=equalities,
        b_eq=np.asarray(eq_rhs),
        bounds=[(None, None)] * (allocation_count + 1),
        method="highs-ipm",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return float(result.x[slack]), len(games), len(edges)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--point-count", type=int, default=8)
    parser.add_argument("--total", type=int, default=10)
    parser.add_argument("--bound", type=int, default=6)
    parser.add_argument("--minimum-nodes", type=int, default=3)
    parser.add_argument("--maximum-components", type=int, default=5)
    parser.add_argument("--bppv-pool", action="store_true")
    parser.add_argument("--bppv-anchor", action="store_true")
    parser.add_argument("--housman-clark-pool", action="store_true")
    parser.add_argument("--housman-clark-anchor", action="store_true")
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.bppv_anchor and not args.bppv_pool:
        parser.error("--bppv-anchor requires --bppv-pool")
    if args.housman_clark_anchor and not args.housman_clark_pool:
        parser.error("--housman-clark-anchor requires --housman-clark-pool")
    if args.bppv_pool and args.housman_clark_pool:
        parser.error("choose only one structured point pool")
    rng = random.Random(args.seed)
    bppv_seed = ((3, 1, 1, 2, 2), (2, 3, 2, 1, 1))
    bppv_pool = sorted(
        {
            tuple(point[index] for index in permutation)
            for point in bppv_seed
            for permutation in itertools.permutations(range(5))
        }
    )
    housman_clark_seed = (
        (0, 0, 4, 4, 0),
        (0, 1, 4, 3, 0),
        (1, 0, 3, 4, 0),
        (3, 4, 1, 0, 0),
        (4, 3, 0, 1, 0),
        (4, 4, 0, 0, 0),
    )
    housman_clark_pool = sorted(
        {
            tuple(point[index] for index in permutation)
            for point in housman_clark_seed
            for permutation in itertools.permutations(range(5))
        }
    )
    best = None
    tested_components = 0
    for trial in range(args.trials):
        if args.bppv_pool:
            if args.n != 5:
                raise ValueError("the BPPV pool is five-player")
            anchors = list(bppv_seed) if args.bppv_anchor else []
            choices = [point for point in bppv_pool if point not in anchors]
            points = anchors + rng.sample(choices, args.point_count - len(anchors))
        elif args.housman_clark_pool:
            if args.n != 5:
                raise ValueError("the Housman-Clark dummy-player pool is five-player")
            anchors = list(housman_clark_seed) if args.housman_clark_anchor else []
            choices = [
                point for point in housman_clark_pool if point not in anchors
            ]
            points = anchors + rng.sample(choices, args.point_count - len(anchors))
        else:
            points = random_points(
                rng, args.n, args.point_count, args.total, args.bound
            )
        games_array, edges = build(points, induced_edges=True)
        components = connected_components(len(games_array), edges)
        for nodes, local_edges in components[: args.maximum_components]:
            if len(nodes) < args.minimum_nodes:
                continue
            games = [tuple(int(value) for value in games_array[node]) for node in nodes]
            margin, orbit_games, orbit_edges = equivariant_margin(
                games, local_edges, args.n
            )
            tested_components += 1
            if best is None or margin < best[0]:
                best = (
                    margin,
                    points,
                    games,
                    local_edges,
                    orbit_games,
                    orbit_edges,
                    trial,
                )
                print(
                    json.dumps(
                        {
                            "trial": trial,
                            "tested_components": tested_components,
                            "base_games": len(games),
                            "base_edges": len(local_edges),
                            "orbit_games": orbit_games,
                            "orbit_edges": orbit_edges,
                            "margin": margin,
                        }
                    ),
                    flush=True,
                )
            if margin < -1e-8:
                break
        if best is not None and best[0] < -1e-8:
            break
    payload = {
        "status": (
            "equivariant_obstruction_found"
            if best is not None and best[0] < -1e-8
            else "no_equivariant_obstruction_found"
        ),
        "n": args.n,
        "configuration": vars(args) | {"output": str(args.output)},
        "tested_components": tested_components,
        "best": None
        if best is None
        else {
            "trial": best[6],
            "margin": best[0],
            "points": [list(point) for point in best[1]],
            "orbit_game_count": best[4],
            "orbit_edge_count": best[5],
            "family": {
                "games": [list(game) for game in best[2]],
                "edges": [list(edge) for edge in best[3]],
            },
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"status": payload["status"], "tested": tested_components, "margin": None if best is None else best[0]}))
    return 1 if payload["status"] == "equivariant_obstruction_found" else 0


if __name__ == "__main__":
    raise SystemExit(main())
