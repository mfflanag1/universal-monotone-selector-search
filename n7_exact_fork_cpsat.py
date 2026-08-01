#!/usr/bin/env python3
"""Search for a three-game exact Housman--Clark fork."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ortools.sat.python import cp_model


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=7)
    parser.add_argument("--game-denominator", type=int, default=3)
    parser.add_argument("--witness-denominator", type=int, default=4)
    parser.add_argument("--time-limit", type=float, default=1800.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.n < 4:
        raise ValueError("at least four players are required")

    n = args.n
    grand = (1 << n) - 1
    game_scale = args.game_denominator
    witness_scale = args.witness_denominator
    payoff_total = game_scale * witness_scale
    model = cp_model.CpModel()

    games = [
        [model.new_int_var(0, game_scale, f"v_{node}_{coalition}") for coalition in range(grand + 1)]
        for node in range(3)
    ]
    for node in range(3):
        model.add(games[node][0] == 0)
        model.add(games[node][grand] == game_scale)
        for coalition in range(grand + 1):
            for player in range(n):
                if not coalition >> player & 1:
                    model.add(games[node][coalition] <= games[node][coalition | (1 << player)])

    changed_a = {}
    changed_b = {}
    for coalition in range(1, grand):
        if coalition & 1:
            changed_a[coalition] = model.new_bool_var(f"changed_a_{coalition}")
        if coalition & 2:
            changed_b[coalition] = model.new_bool_var(f"changed_b_{coalition}")
    model.add(sum(changed_a.values()) == 1)
    model.add(sum(changed_b.values()) == 1)

    for coalition in range(grand + 1):
        if coalition in changed_a:
            flag = changed_a[coalition]
            model.add(games[1][coalition] - games[0][coalition] >= flag)
            model.add(games[1][coalition] - games[0][coalition] <= game_scale * flag)
        else:
            model.add(games[1][coalition] == games[0][coalition])
        if coalition in changed_b:
            flag = changed_b[coalition]
            model.add(games[2][coalition] - games[0][coalition] >= flag)
            model.add(games[2][coalition] - games[0][coalition] <= game_scale * flag)
        else:
            model.add(games[2][coalition] == games[0][coalition])

    # A game is exact iff every coalition has a core allocation tight on it.
    for node in range(3):
        for tight_coalition in range(1, grand):
            witness = [
                model.new_int_var(0, payoff_total, f"x_{node}_{tight_coalition}_{player}")
                for player in range(n)
            ]
            model.add(sum(witness) == payoff_total)
            for coalition in range(1, grand):
                coalition_sum = sum(
                    witness[player] for player in range(n) if coalition >> player & 1
                )
                model.add(coalition_sum >= witness_scale * games[node][coalition])
                if coalition == tight_coalition:
                    model.add(coalition_sum == witness_scale * games[node][coalition])

    # Exactness turns the three terminal lower expectations into game values:
    # min_U(x0+x1) - max_V(x0) - max_W(x1).
    obstruction = (
        games[0][3]
        + games[1][grand ^ 1]
        + games[2][grand ^ 2]
        - 2 * game_scale
    )
    model.add(obstruction >= 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.time_limit
    solver.parameters.num_search_workers = args.workers
    solver.parameters.log_search_progress = True
    status = solver.solve(model)
    payload = {
        "status": solver.status_name(status),
        "n": n,
        "game_denominator": game_scale,
        "witness_denominator": witness_scale,
        "wall_time": solver.wall_time,
    }
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        payload.update(
            {
                "obstruction_numerator": solver.value(obstruction),
                "changed_coalitions": [
                    next(coalition for coalition, flag in changed_a.items() if solver.value(flag)),
                    next(coalition for coalition, flag in changed_b.items() if solver.value(flag)),
                ],
                "games": [
                    [solver.value(value) for value in game]
                    for game in games
                ],
            }
        )
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key != "games"}, indent=2))
    return 1 if payload.get("obstruction_numerator", 0) > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
