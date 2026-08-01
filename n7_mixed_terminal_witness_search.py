#!/usr/bin/env python3
"""Mixed-terminal search with exactness enforced by tight core witnesses."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from n5_mixed_terminal_search import (
    Topology,
    lower_expectation_dual,
    random_topology,
    sparse_matrix,
)


class WitnessExactCone:
    def __init__(self, n: int, topology: Topology) -> None:
        self.n = n
        self.grand = (1 << n) - 1
        self.coalition_count = self.grand + 1
        self.proper = tuple(range(1, self.grand))
        self.node_count = topology.node_count
        self.game_variable_count = self.node_count * self.coalition_count
        self.witness_offset = self.game_variable_count
        self.variable_count = self.witness_offset + (
            self.node_count * len(self.proper) * n
        )
        inequalities: list[dict[int, float]] = []
        inequality_rhs: list[float] = []
        equalities: list[dict[int, float]] = []
        equality_rhs: list[float] = []

        for node in range(self.node_count):
            equalities.append({self.game_column(node, 0): 1.0})
            equality_rhs.append(0.0)
            equalities.append({self.game_column(node, self.grand): 1.0})
            equality_rhs.append(1.0)
            for coalition in range(self.coalition_count):
                for player in range(n):
                    if coalition >> player & 1:
                        continue
                    successor = coalition | (1 << player)
                    inequalities.append(
                        {
                            self.game_column(node, coalition): 1.0,
                            self.game_column(node, successor): -1.0,
                        }
                    )
                    inequality_rhs.append(0.0)
            for tight in self.proper:
                efficiency = {
                    self.witness_column(node, tight, player): 1.0
                    for player in range(n)
                }
                efficiency[self.game_column(node, self.grand)] = -1.0
                equalities.append(efficiency)
                equality_rhs.append(0.0)
                tightness = {
                    self.game_column(node, tight): 1.0,
                    **{
                        self.witness_column(node, tight, player): -1.0
                        for player in range(n)
                        if tight >> player & 1
                    },
                }
                equalities.append(tightness)
                equality_rhs.append(0.0)
                for coalition in self.proper:
                    row = {self.game_column(node, coalition): 1.0}
                    for player in range(n):
                        if coalition >> player & 1:
                            column = self.witness_column(node, tight, player)
                            row[column] = row.get(column, 0.0) - 1.0
                    inequalities.append(row)
                    inequality_rhs.append(0.0)

        for lower, upper, protected in topology.arcs:
            for coalition in range(self.coalition_count):
                row = {
                    self.game_column(lower, coalition): 1.0,
                    self.game_column(upper, coalition): -1.0,
                }
                if coalition & protected != protected:
                    equalities.append(row)
                    equality_rhs.append(0.0)
                else:
                    inequalities.append(row)
                    inequality_rhs.append(0.0)

        self.a_ub = sparse_matrix(inequalities, self.variable_count)
        self.b_ub = np.asarray(inequality_rhs)
        self.a_eq = sparse_matrix(equalities, self.variable_count)
        self.b_eq = np.asarray(equality_rhs)

    def game_column(self, node: int, coalition: int) -> int:
        return node * self.coalition_count + coalition

    def witness_column(self, node: int, tight: int, player: int) -> int:
        return self.witness_offset + (
            (node * len(self.proper) + tight - 1) * self.n + player
        )

    def maximize(self, game_objective: np.ndarray) -> np.ndarray | None:
        objective = np.zeros(self.variable_count)
        objective[: self.game_variable_count] = -game_objective
        result = linprog(
            objective,
            A_ub=self.a_ub,
            b_ub=self.b_ub,
            A_eq=self.a_eq,
            b_eq=self.b_eq,
            bounds=[(0.0, 1.0)] * self.variable_count,
            method="highs",
        )
        return result.x if result.success else None

    def common_point_seed(self, points: list[tuple[float, ...]]) -> np.ndarray:
        game = np.asarray(
            [
                min(
                    sum(point[player] for player in range(self.n) if coalition >> player & 1)
                    for point in points
                )
                for coalition in range(self.coalition_count)
            ]
        )
        seed = np.zeros(self.variable_count)
        for node in range(self.node_count):
            seed[
                node * self.coalition_count : (node + 1) * self.coalition_count
            ] = game
            for tight in self.proper:
                witness = min(
                    points,
                    key=lambda point: sum(
                        point[player]
                        for player in range(self.n)
                        if tight >> player & 1
                    ),
                )
                for player, value in enumerate(witness):
                    seed[self.witness_column(node, tight, player)] = value
        return seed


def optimize(
    topology: Topology,
    n: int,
    rng: random.Random,
    starts: int,
    iterations: int,
) -> dict | None:
    cone = WitnessExactCone(n, topology)
    best = None
    initial_points = []
    if n >= 5:
        first = (3, 1, 1, 2, 2)
        second = (2, 3, 2, 1, 1)
        for permutation in itertools.islice(itertools.permutations(range(5)), starts):
            points = [
                tuple(seed[permutation[player]] / 9.0 for player in range(5))
                + (0.0,) * (n - 5)
                for seed in (first, second)
            ]
            initial_points.append(cone.common_point_seed(points))
    for _ in range(starts):
        point = cone.maximize(
            np.asarray(
                [
                    rng.uniform(-1.0, 1.0)
                    for _ in range(cone.game_variable_count)
                ]
            )
        )
        if point is None:
            continue
        initial_points.append(point)
    for point in initial_points:
        for _ in range(iterations):
            objective = np.zeros(cone.game_variable_count)
            valid = True
            for node, divergence in enumerate(topology.divergence):
                game = point[
                    node * cone.coalition_count :
                    (node + 1) * cone.coalition_count
                ]
                dual = lower_expectation_dual(game, divergence, n)
                if dual is None:
                    valid = False
                    break
                _, weights, _ = dual
                for index, coalition in enumerate(cone.proper):
                    objective[cone.game_column(node, coalition)] += weights[index]
            if not valid:
                break
            updated = cone.maximize(objective)
            if updated is None:
                break
            point = updated
        values = []
        for node, divergence in enumerate(topology.divergence):
            game = point[
                node * cone.coalition_count :
                (node + 1) * cone.coalition_count
            ]
            dual = lower_expectation_dual(game, divergence, n)
            if dual is None:
                values = []
                break
            values.append(dual[0])
        if not values:
            continue
        candidate = {
            "objective": sum(values),
            "node_lower_expectations": values,
            "games": [
                [
                    float(value)
                    for value in point[
                        node * cone.coalition_count :
                        (node + 1) * cone.coalition_count
                    ]
                ]
                for node in range(topology.node_count)
            ],
        }
        if best is None or candidate["objective"] > best["objective"]:
            best = candidate
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=7)
    parser.add_argument("--nodes", type=int, default=5)
    parser.add_argument("--samples", type=int, default=20)
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    rows = []
    tested = 0
    while tested < args.samples:
        topology = random_topology(args.n, args.nodes, rng)
        if topology is None:
            continue
        tested += 1
        result = optimize(
            topology, args.n, rng, args.starts, args.iterations
        )
        if result is None:
            continue
        row = {
            "arcs": [list(arc) for arc in topology.arcs],
            "divergence": [list(vector) for vector in topology.divergence],
            **result,
        }
        rows.append(row)
        rows.sort(key=lambda item: item["objective"], reverse=True)
        del rows[10:]
        print(
            json.dumps(
                {"tested": tested, "objective": result["objective"]}
            ),
            flush=True,
        )
        if result["objective"] > 1e-8:
            break
    payload = {
        "status": (
            "positive_terminal_obstruction_found"
            if rows and rows[0]["objective"] > 1e-8
            else "no_positive_terminal_obstruction_found"
        ),
        "n": args.n,
        "node_count": args.nodes,
        "tested_topologies": tested,
        "best": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in payload if key != "best"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
