#!/usr/bin/env python3
"""Globally optimize all unit-indicator k-by-k terminal routings at n=6."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from n6_mixed_terminal_search import GRAND, N, load_facets


class Model:
    def __init__(self, k: int) -> None:
        self.k = k
        self.proper = tuple(range(1, GRAND + 1))
        self.node_count = 2 * k
        self.game_count = self.node_count * (GRAND + 1)
        self.z_source_offset = self.game_count
        self.z_sink_offset = self.z_source_offset + k * len(self.proper)
        self.route_offset = self.z_sink_offset + k * len(self.proper)
        self.source_value_offset = self.route_offset + k * k * N
        self.sink_value_offset = self.source_value_offset + k
        self.variable_count = self.sink_value_offset + k
        self.rows: list[dict[int, float]] = []
        self.lower: list[float] = []
        self.upper: list[float] = []

    def game(self, node: int, coalition: int) -> int:
        return node * (GRAND + 1) + coalition

    def z_source(self, source: int, coalition: int) -> int:
        return self.z_source_offset + source * len(self.proper) + coalition - 1

    def z_sink(self, sink: int, coalition: int) -> int:
        return self.z_sink_offset + sink * len(self.proper) + coalition - 1

    def route(self, source: int, sink: int, player: int) -> int:
        return self.route_offset + (source * self.k + sink) * N + player

    def add(self, row: dict[int, float], lower=-np.inf, upper=np.inf) -> None:
        self.rows.append({column: value for column, value in row.items() if value})
        self.lower.append(lower)
        self.upper.append(upper)

    def matrix(self):
        rr, cc, vv = [], [], []
        for row_index, row in enumerate(self.rows):
            for column, value in row.items():
                rr.append(row_index)
                cc.append(column)
                vv.append(value)
        return coo_matrix(
            (vv, (rr, cc)), shape=(len(self.rows), self.variable_count)
        ).tocsr()


def build(k: int) -> tuple[Model, np.ndarray, np.ndarray, Bounds, LinearConstraint]:
    model = Model(k)
    facets = load_facets()
    for node in range(model.node_count):
        model.add({model.game(node, 0): 1.0}, 0.0, 0.0)
        model.add({model.game(node, GRAND): 1.0}, 1.0, 1.0)
        for facet in facets:
            model.add(
                {
                    model.game(node, coalition): float(coefficient)
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                },
                0.0,
            )
        for coalition in range(GRAND + 1):
            for player in range(N):
                if coalition >> player & 1:
                    continue
                successor = coalition | (1 << player)
                model.add(
                    {
                        model.game(node, successor): 1.0,
                        model.game(node, coalition): -1.0,
                    },
                    0.0,
                )

    for source in range(k):
        model.add(
            {model.z_source(source, coalition): 1.0 for coalition in model.proper},
            1.0,
            1.0,
        )
    for sink in range(k):
        model.add(
            {model.z_sink(sink, coalition): 1.0 for coalition in model.proper},
            1.0,
            1.0,
        )

    for source in range(k):
        for player in range(N):
            row = {model.route(source, sink, player): 1.0 for sink in range(k)}
            for coalition in model.proper:
                if coalition >> player & 1:
                    row[model.z_source(source, coalition)] = -1.0
            model.add(row, 0.0, 0.0)
    for sink in range(k):
        for player in range(N):
            row = {model.route(source, sink, player): 1.0 for source in range(k)}
            for coalition in model.proper:
                if coalition >> player & 1:
                    row[model.z_sink(sink, coalition)] = -1.0
            model.add(row, 0.0, 0.0)

    for source in range(k):
        selected = model.source_value_offset + source
        for coalition in model.proper:
            binary = model.z_source(source, coalition)
            game = model.game(source, coalition)
            model.add({selected: 1.0, game: -1.0, binary: 1.0}, upper=1.0)
            model.add({selected: 1.0, game: -1.0, binary: -1.0}, lower=-1.0)
    for sink in range(k):
        selected = model.sink_value_offset + sink
        node = k + sink
        for coalition in model.proper:
            binary = model.z_sink(sink, coalition)
            game = model.game(node, GRAND ^ coalition)
            model.add({selected: 1.0, game: -1.0, binary: 1.0}, upper=1.0)
            model.add({selected: 1.0, game: -1.0, binary: -1.0}, lower=-1.0)

    for source in range(k):
        for sink in range(k):
            upper_node = k + sink
            for player in range(N):
                binary = model.route(source, sink, player)
                for coalition in range(GRAND + 1):
                    difference = {
                        model.game(source, coalition): 1.0,
                        model.game(upper_node, coalition): -1.0,
                        binary: 1.0,
                    }
                    model.add(difference, upper=1.0)
                    if not coalition >> player & 1:
                        reverse = {
                            model.game(source, coalition): 1.0,
                            model.game(upper_node, coalition): -1.0,
                            binary: -1.0,
                        }
                        model.add(reverse, lower=-1.0)

    objective = np.zeros(model.variable_count)
    objective[model.source_value_offset : model.source_value_offset + k] = -1.0
    objective[model.sink_value_offset : model.sink_value_offset + k] = -1.0
    integrality = np.zeros(model.variable_count, dtype=np.uint8)
    integrality[model.z_source_offset : model.source_value_offset] = 1
    lower = np.zeros(model.variable_count)
    upper = np.ones(model.variable_count)
    return (
        model,
        objective,
        integrality,
        Bounds(lower, upper),
        LinearConstraint(model.matrix(), np.asarray(model.lower), np.asarray(model.upper)),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=3)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    model, objective, integrality, bounds, constraints = build(args.k)
    result = milp(
        objective,
        integrality=integrality,
        bounds=bounds,
        constraints=constraints,
        options={"time_limit": args.time_limit, "mip_rel_gap": 0.0},
    )
    selected_sum = None if result.fun is None else -float(result.fun)
    terminal_objective = None if selected_sum is None else selected_sum - args.k
    payload = {
        "status": result.message,
        "success": bool(result.success),
        "mip_gap": getattr(result, "mip_gap", None),
        "mip_node_count": getattr(result, "mip_node_count", None),
        "selected_worth_sum": selected_sum,
        "terminal_objective": terminal_objective,
        "variable_count": model.variable_count,
        "constraint_count": len(model.rows),
    }
    if result.x is not None:
        sources = []
        sinks = []
        supports = [[0] * args.k for _ in range(args.k)]
        for source in range(args.k):
            sources.append(
                max(model.proper, key=lambda coalition: result.x[model.z_source(source, coalition)])
            )
        for sink in range(args.k):
            sinks.append(
                max(model.proper, key=lambda coalition: result.x[model.z_sink(sink, coalition)])
            )
        for source in range(args.k):
            for sink in range(args.k):
                for player in range(N):
                    if result.x[model.route(source, sink, player)] > 0.5:
                        supports[source][sink] |= 1 << player
        payload.update(
            {
                "sources": sources,
                "sinks": sinks,
                "supports": supports,
                "games": [
                    [
                        float(value)
                        for value in result.x[
                            node * (GRAND + 1) : (node + 1) * (GRAND + 1)
                        ]
                    ]
                    for node in range(model.node_count)
                ],
            }
        )
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key != "games"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
