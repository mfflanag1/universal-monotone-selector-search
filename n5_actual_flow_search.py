#!/usr/bin/env python3
"""Search explicit protected flows on legal exact five-player families."""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_facet_search import ExactConeModel, random_objective, solve_model
from n5_actual_bowtie_margin_search import games_from_point, solve_minimum
from n5_actual_bowtie_milp import FlowTopology, LegalExactCone
from n5_actual_bowtie_search import matrices
from n5_mixed_terminal_search import lower_expectation_dual
from n5_random_exact_hourglass import hourglass


N = 5
GRAND = (1 << N) - 1


def random_flow(
    edges: list[tuple[int, int, int]], rng: random.Random
) -> tuple[list[list[int]], list[list[int]]]:
    node_count = 1 + max(max(lower, upper) for lower, upper, _ in edges)
    divergence = [[0] * N for _ in range(node_count)]
    flows = []
    for lower, upper, protected in edges:
        flow = [
            rng.randint(1, 5) if protected >> player & 1 else 0
            for player in range(N)
        ]
        flows.append(flow)
        for player, weight in enumerate(flow):
            divergence[lower][player] += weight
            divergence[upper][player] -= weight
    return flows, divergence


def optimize_flow(
    edges: list[tuple[int, int, int]],
    bumps: list[float],
    divergence: list[list[int]],
    starts: int,
    iterations: int,
    rng: random.Random,
    tangent_player: int | None = None,
) -> dict[str, Any] | None:
    if tangent_player is None:
        model = ExactConeModel(edges, bumps)
        variable_count = model.variable_count
        game_column = model.game_index

        def initial_point():
            return solve_model(model, random_objective(model, rng))

        def optimize_objective(objective):
            return solve_model(model, objective)

        def extract_games(point):
            return point

    else:
        topology = FlowTopology(
            len(divergence), tuple(edges), tuple(), tuple()
        )
        model = LegalExactCone(
            topology, min(bumps), tangent_player=tangent_player
        )
        matrix_data = matrices(model)
        variable_count = model.game_variables
        game_column = model.column

        def initial_point():
            objective = np.asarray(
                [rng.uniform(-1.0, 1.0) for _ in range(variable_count)]
            )
            return solve_minimum(model, matrix_data, objective)

        def optimize_objective(objective):
            return solve_minimum(model, matrix_data, objective)

        def extract_games(point):
            return games_from_point(model, point)

    best = None
    for _ in range(starts):
        point = initial_point()
        if point is None:
            continue
        games = extract_games(point)
        for _ in range(iterations):
            objective = np.zeros(variable_count)
            valid = True
            for node, vector in enumerate(divergence):
                dual = lower_expectation_dual(
                    np.asarray(games[node]), tuple(vector), N
                )
                if dual is None:
                    valid = False
                    break
                _value, weights, beta = dual
                for coalition, weight in enumerate(weights, start=1):
                    if weight:
                        objective[game_column(node, coalition)] -= weight
                objective[game_column(node, GRAND)] -= beta
            if not valid:
                break
            point = optimize_objective(objective)
            if point is None:
                break
            games = extract_games(point)
        values = []
        for game, vector in zip(games, divergence, strict=True):
            dual = lower_expectation_dual(
                np.asarray(game), tuple(vector), N
            )
            if dual is None:
                values = []
                break
            values.append(dual[0])
        if not values:
            continue
        total = sum(values)
        if best is None or total > float(best["objective"]):
            best = {
                "objective": total,
                "node_lower_expectations": values,
                "games": games,
            }
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=100)
    parser.add_argument("--starts", type=int, default=20)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--bump", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--tangent-player", type=int, choices=range(-1, N))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pool = [mask for mask in range(1, 31) if 2 <= mask.bit_count() <= 4]
    best = []
    for case in range(args.cases):
        incoming = rng.sample(pool, rng.randint(2, 6))
        outgoing = rng.sample(pool, rng.randint(2, 6))
        edges = hourglass(incoming, outgoing)
        bumps = [args.bump] * len(edges)
        flows, divergence = random_flow(edges, rng)
        optimum = optimize_flow(
            edges,
            bumps,
            divergence,
            args.starts,
            args.iterations,
            rng,
            args.tangent_player,
        )
        if optimum is None:
            continue
        row = {
            "case": case,
            "incoming": incoming,
            "outgoing": outgoing,
            "edges": [list(edge) for edge in edges],
            "flows": flows,
            "divergence": divergence,
            **optimum,
        }
        best.append(row)
        best.sort(key=lambda item: float(item["objective"]), reverse=True)
        del best[20:]
        print(
            json.dumps({"case": case, "objective": row["objective"]}),
            flush=True,
        )
        if float(row["objective"]) > 1e-8:
            break
    payload = {
        "status": (
            "positive_flow_found"
            if best and float(best[0]["objective"]) > 1e-8
            else "no_positive_flow_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "best": best,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"status": payload["status"], "best": best[:1]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
