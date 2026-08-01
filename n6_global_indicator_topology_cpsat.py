#!/usr/bin/env python3
"""Exact integer search over all n=6 unit-indicator terminal routings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ortools.sat.python import cp_model

from n6_mixed_terminal_search import GRAND, N, load_facets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=3)
    parser.add_argument("--denominator", type=int, default=4)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--source0-size", type=int)
    parser.add_argument("--positive-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    model = cp_model.CpModel()
    coalitions = tuple(range(1, GRAND + 1))
    node_count = 2 * args.k
    games = [
        [
            model.new_int_var(0, args.denominator, f"v_{node}_{coalition}")
            for coalition in range(GRAND + 1)
        ]
        for node in range(node_count)
    ]
    facets = load_facets()
    for node in range(node_count):
        model.add(games[node][0] == 0)
        model.add(games[node][GRAND] == args.denominator)
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
                if coalition >> player & 1:
                    continue
                model.add(games[node][coalition] <= games[node][coalition | (1 << player)])

    source_choice = [
        [model.new_bool_var(f"zs_{source}_{coalition}") for coalition in coalitions]
        for source in range(args.k)
    ]
    sink_choice = [
        [model.new_bool_var(f"zt_{sink}_{coalition}") for coalition in coalitions]
        for sink in range(args.k)
    ]
    for choices in (*source_choice, *sink_choice):
        model.add_exactly_one(choices)
    if args.source0_size is None:
        for source in range(args.k - 1):
            model.add(
                sum(coalition * source_choice[source][coalition - 1] for coalition in coalitions)
                <= sum(coalition * source_choice[source + 1][coalition - 1] for coalition in coalitions)
            )
    for sink in range(args.k - 1):
        model.add(
            sum(coalition * sink_choice[sink][coalition - 1] for coalition in coalitions)
            <= sum(coalition * sink_choice[sink + 1][coalition - 1] for coalition in coalitions)
        )
    if args.source0_size is not None:
        if args.source0_size not in range(1, N + 1):
            raise ValueError("source0 size must be between one and six")
        canonical = (1 << args.source0_size) - 1
        model.add(source_choice[0][canonical - 1] == 1)

    route = [
        [
            [model.new_bool_var(f"y_{source}_{sink}_{player}") for player in range(N)]
            for sink in range(args.k)
        ]
        for source in range(args.k)
    ]
    for source in range(args.k):
        for player in range(N):
            model.add(
                sum(route[source][sink][player] for sink in range(args.k))
                == sum(
                    source_choice[source][coalition - 1]
                    for coalition in coalitions
                    if coalition >> player & 1
                )
            )
    for sink in range(args.k):
        for player in range(N):
            model.add(
                sum(route[source][sink][player] for source in range(args.k))
                == sum(
                    sink_choice[sink][coalition - 1]
                    for coalition in coalitions
                    if coalition >> player & 1
                )
            )

    source_values = [
        model.new_int_var(0, args.denominator, f"source_value_{source}")
        for source in range(args.k)
    ]
    sink_values = [
        model.new_int_var(0, args.denominator, f"sink_value_{sink}")
        for sink in range(args.k)
    ]
    for source in range(args.k):
        for coalition in coalitions:
            model.add(source_values[source] == games[source][coalition]).only_enforce_if(
                source_choice[source][coalition - 1]
            )
    for sink in range(args.k):
        node = args.k + sink
        for coalition in coalitions:
            model.add(sink_values[sink] == games[node][GRAND ^ coalition]).only_enforce_if(
                sink_choice[sink][coalition - 1]
            )

    for source in range(args.k):
        for sink in range(args.k):
            node = args.k + sink
            for player in range(N):
                active = route[source][sink][player]
                for coalition in range(GRAND + 1):
                    model.add(games[source][coalition] <= games[node][coalition]).only_enforce_if(active)
                    if not coalition >> player & 1:
                        model.add(games[source][coalition] == games[node][coalition]).only_enforce_if(active)

    selected_expression = sum(source_values) + sum(sink_values)
    if args.positive_only:
        model.add(selected_expression >= args.k * args.denominator + 1)
    else:
        model.maximize(selected_expression)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.time_limit
    solver.parameters.num_search_workers = args.workers
    solver.parameters.log_search_progress = True
    status = solver.solve(model)
    selected_sum = None
    terminal_objective = None
    payload = {
        "status": solver.status_name(status),
        "n": N,
        "k": args.k,
        "denominator": args.denominator,
        "facet_count": len(facets),
        "positive_only": args.positive_only,
        "wall_time": solver.wall_time,
        "best_bound": solver.best_objective_bound,
    }
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        selected_sum = sum(solver.value(value) for value in (*source_values, *sink_values))
        terminal_objective = selected_sum - args.k * args.denominator
        sources = []
        sinks = []
        supports = [[0] * args.k for _ in range(args.k)]
        for source in range(args.k):
            sources.append(
                next(
                    coalition
                    for coalition in coalitions
                    if solver.value(source_choice[source][coalition - 1])
                )
            )
        for sink in range(args.k):
            sinks.append(
                next(
                    coalition
                    for coalition in coalitions
                    if solver.value(sink_choice[sink][coalition - 1])
                )
            )
        for source in range(args.k):
            for sink in range(args.k):
                for player in range(N):
                    if solver.value(route[source][sink][player]):
                        supports[source][sink] |= 1 << player
        payload.update(
            {
                "selected_worth_sum": selected_sum,
                "terminal_objective_numerator": terminal_objective,
                "sources": sources,
                "sinks": sinks,
                "supports": supports,
                "games": [
                    [solver.value(games[node][coalition]) for coalition in range(GRAND + 1)]
                    for node in range(node_count)
                ],
            }
        )
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key != "games"}, indent=2))
    return 1 if terminal_objective is not None and terminal_objective > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
