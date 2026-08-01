#!/usr/bin/env python3
"""Search the complete n=6 exact cone for a mixed protected-flow obstruction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ortools.sat.python import cp_model

from n6_mixed_terminal_search import load_facets as load_n6_facets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(5, 6), default=6)
    parser.add_argument("--nodes", type=int, default=4)
    parser.add_argument("--game-denominator", type=int, default=4)
    parser.add_argument("--grand-maximum", type=int, default=0)
    parser.add_argument("--payoff-denominator", type=int, default=4)
    parser.add_argument("--flow-bound", type=int, default=3)
    parser.add_argument("--normal-bound", type=int, default=6)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--hourglass", action="store_true")
    parser.add_argument("--mixed-fork", action="store_true")
    parser.add_argument("--complementary-mixed-fork", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.nodes < 2:
        raise ValueError("at least two terminal nodes are required")
    if args.hourglass and args.nodes != 5:
        raise ValueError("the hourglass topology requires five nodes")
    if (args.mixed_fork or args.complementary_mixed_fork) and args.nodes != 3:
        raise ValueError("the mixed-fork topology requires three nodes")
    if sum((args.hourglass, args.mixed_fork, args.complementary_mixed_fork)) > 1:
        raise ValueError("choose at most one fixed topology")
    if (args.mixed_fork or args.complementary_mixed_fork) and (args.flow_bound < 2 or args.n < 3):
        raise ValueError("the mixed fork fixes two player flows to 1 and 2")
    N = args.n
    GRAND = (1 << N) - 1
    game_maximum = args.grand_maximum or args.game_denominator
    if game_maximum < args.game_denominator:
        raise ValueError("grand maximum cannot be below the game denominator")
    if N == 5:
        from n5_facet_search import load_facets as load_n5_facets

        facets, _types = load_n5_facets()
    else:
        facets = load_n6_facets()

    model = cp_model.CpModel()
    proper = range(1, GRAND)
    games = [
        [
            model.new_int_var(0, game_maximum, f"v_{node}_{coalition}")
            for coalition in range(GRAND + 1)
        ]
        for node in range(args.nodes)
    ]
    for node in range(args.nodes):
        model.add(games[node][0] == 0)
        if game_maximum == args.game_denominator:
            model.add(games[node][GRAND] == args.game_denominator)
        else:
            model.add(games[node][GRAND] >= args.game_denominator)
        for facet in facets:
            model.add(
                sum(
                    coefficient * games[node][coalition]
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                )
                >= 0
            )
        for coalition in range(GRAND + 1):
            for player in range(N):
                if not coalition >> player & 1:
                    model.add(games[node][coalition] <= games[node][coalition | (1 << player)])

    payoff_limit = args.payoff_denominator * game_maximum
    payoffs = [
        [model.new_int_var(0, payoff_limit, f"x_{node}_{player}") for player in range(N)]
        for node in range(args.nodes)
    ]
    for node in range(args.nodes):
        model.add(
            sum(payoffs[node])
            == args.payoff_denominator * games[node][GRAND]
        )
        for coalition in proper:
            model.add(
                sum(payoffs[node][player] for player in range(N) if coalition >> player & 1)
                >= args.payoff_denominator * games[node][coalition]
            )

    flows = {}
    active_routes = {}
    for lower in range(args.nodes):
        for upper in range(args.nodes):
            if lower == upper:
                continue
            for player in range(N):
                key = (lower, upper, player)
                flow = model.new_int_var(0, args.flow_bound, f"flow_{lower}_{upper}_{player}")
                active = model.new_bool_var(f"active_{lower}_{upper}_{player}")
                flows[key] = flow
                active_routes[key] = active
                model.add(flow >= active)
                model.add(flow <= args.flow_bound * active)
                allowed_pairs = (
                    {(0, 2), (1, 2), (2, 3), (2, 4)}
                    if args.hourglass
                    else {(0, 1), (0, 2)}
                    if args.mixed_fork or args.complementary_mixed_fork
                    else None
                )
                if allowed_pairs is not None and (lower, upper) not in allowed_pairs:
                    model.add(active == 0)
                for coalition in range(GRAND + 1):
                    model.add(games[lower][coalition] <= games[upper][coalition]).only_enforce_if(active)
                    if not coalition >> player & 1:
                        model.add(games[lower][coalition] == games[upper][coalition]).only_enforce_if(active)

    if args.hourglass:
        for lower, upper in ((0, 2), (1, 2), (2, 3), (2, 4)):
            model.add(sum(active_routes[lower, upper, player] for player in range(N)) >= 1)
    if args.mixed_fork or args.complementary_mixed_fork:
        for upper in (1, 2):
            model.add(sum(active_routes[0, upper, player] for player in range(N)) >= 1)
        model.add(flows[0, 1, 0] == 1)
        model.add(flows[0, 1, 1] == 2)
    if args.complementary_mixed_fork:
        for player in range(2, N):
            model.add(flows[0, 1, player] == 0)
        model.add(flows[0, 2, 0] == 1)
        model.add(flows[0, 2, 1] == 0)
        for player in range(2, N):
            model.add(flows[0, 2, player] == 2)

    divergence_bound = (args.nodes - 1) * args.flow_bound
    divergences = [
        [
            model.new_int_var(-divergence_bound, divergence_bound, f"d_{node}_{player}")
            for player in range(N)
        ]
        for node in range(args.nodes)
    ]
    for node in range(args.nodes):
        for player in range(N):
            model.add(
                divergences[node][player]
                == sum(
                    flows[node, upper, player]
                    for upper in range(args.nodes)
                    if upper != node
                )
                - sum(
                    flows[lower, node, player]
                    for lower in range(args.nodes)
                    if lower != node
                )
            )
    if args.hourglass:
        central_positive = [model.new_bool_var(f"central_positive_{player}") for player in range(N)]
        central_negative = [model.new_bool_var(f"central_negative_{player}") for player in range(N)]
        for player in range(N):
            model.add(divergences[2][player] >= 1).only_enforce_if(central_positive[player])
            model.add(divergences[2][player] <= 0).only_enforce_if(central_positive[player].Not())
            model.add(divergences[2][player] <= -1).only_enforce_if(central_negative[player])
            model.add(divergences[2][player] >= 0).only_enforce_if(central_negative[player].Not())
        model.add(sum(central_positive) >= 1)
        model.add(sum(central_negative) >= 1)
    # Every nonzero obstruction can be relabeled to have this positive entry.
    model.add(divergences[0][0] >= 1)

    normal_multipliers = []
    normal_constants = []
    for node in range(args.nodes):
        multipliers = [
            model.new_int_var(0, args.normal_bound, f"lambda_{node}_{coalition}")
            for coalition in proper
        ]
        normal_multipliers.append(multipliers)
        alpha = model.new_int_var(
            -args.normal_bound, args.normal_bound, f"alpha_{node}"
        )
        normal_constants.append(alpha)
        for coalition in proper:
            positive = model.new_bool_var(f"lambda_positive_{node}_{coalition}")
            multiplier = multipliers[coalition - 1]
            model.add(multiplier >= positive)
            model.add(multiplier <= args.normal_bound * positive)
            model.add(
                sum(payoffs[node][player] for player in range(N) if coalition >> player & 1)
                == args.payoff_denominator * games[node][coalition]
            ).only_enforce_if(positive)
        for player in range(N):
            model.add(
                divergences[node][player]
                == alpha
                + sum(
                    multipliers[coalition - 1]
                    for coalition in proper
                    if coalition >> player & 1
                )
            )

    products = []
    for node in range(args.nodes):
        for player in range(N):
            product = model.new_int_var(
                -divergence_bound * payoff_limit,
                divergence_bound * payoff_limit,
                f"objective_product_{node}_{player}",
            )
            model.add_multiplication_equality(
                product, [divergences[node][player], payoffs[node][player]]
            )
            products.append(product)
    model.add(sum(products) >= 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.time_limit
    solver.parameters.num_search_workers = args.workers
    solver.parameters.log_search_progress = True
    status = solver.solve(model)
    payload = {
        "status": solver.status_name(status),
        "n": N,
        "nodes": args.nodes,
        "game_denominator": args.game_denominator,
        "grand_maximum": game_maximum,
        "payoff_denominator": args.payoff_denominator,
        "flow_bound": args.flow_bound,
        "normal_bound": args.normal_bound,
        "facet_count": len(facets),
        "hourglass": args.hourglass,
        "mixed_fork": args.mixed_fork,
        "complementary_mixed_fork": args.complementary_mixed_fork,
        "wall_time": solver.wall_time,
    }
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        game_values = [
            [solver.value(value) for value in game]
            for game in games
        ]
        payoff_values = [
            [solver.value(value) for value in payoff]
            for payoff in payoffs
        ]
        divergence_values = [
            [solver.value(value) for value in divergence]
            for divergence in divergences
        ]
        flow_values = [
            [lower, upper, player, solver.value(flow)]
            for (lower, upper, player), flow in flows.items()
            if solver.value(flow)
        ]
        multiplier_values = [
            {
                str(coalition): solver.value(normal_multipliers[node][coalition - 1])
                for coalition in proper
                if solver.value(normal_multipliers[node][coalition - 1])
            }
            for node in range(args.nodes)
        ]
        objective_numerator = sum(
            divergence_values[node][player] * payoff_values[node][player]
            for node in range(args.nodes)
            for player in range(N)
        )
        payload.update(
            {
                "objective_numerator": objective_numerator,
                "games": game_values,
                "payoffs": payoff_values,
                "divergences": divergence_values,
                "flows": flow_values,
                "normal_constants": [solver.value(value) for value in normal_constants],
                "normal_multipliers": multiplier_values,
            }
        )
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key not in {"games", "normal_multipliers"}}, indent=2))
    return 1 if payload.get("objective_numerator", 0) > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
