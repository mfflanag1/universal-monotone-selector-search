#!/usr/bin/env python3
"""Coordinate-ascent search for mixed-divergence terminal obstructions."""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

sys.path.append(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)

from n5_facet_search import load_facets


def sparse_matrix(rows: list[dict[int, float]], columns: int):
    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[float] = []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                row_indices.append(row_index)
                column_indices.append(column)
                values.append(value)
    return coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(rows), columns),
    ).tocsr()


@dataclass(frozen=True)
class Topology:
    node_count: int
    arcs: tuple[tuple[int, int, int], ...]
    divergence: tuple[tuple[int, ...], ...]


class GameCone:
    def __init__(self, n: int, topology: Topology) -> None:
        self.n = n
        self.grand = (1 << n) - 1
        self.coalition_count = self.grand + 1
        self.node_count = topology.node_count
        self.variable_count = self.node_count * self.coalition_count
        self.inequalities: list[dict[int, float]] = []
        self.inequality_rhs: list[float] = []
        self.equalities: list[dict[int, float]] = []
        self.equality_rhs: list[float] = []
        facets, _ = load_facets()
        for node in range(self.node_count):
            self.equalities.append({self.column(node, 0): 1.0})
            self.equality_rhs.append(0.0)
            self.equalities.append({self.column(node, self.grand): 1.0})
            self.equality_rhs.append(1.0)
            for facet in facets:
                self.inequalities.append(
                    {
                        self.column(node, coalition): -float(coefficient)
                        for coalition, coefficient in enumerate(facet)
                        if coefficient
                    }
                )
                self.inequality_rhs.append(0.0)
            for coalition in range(self.coalition_count):
                for player in range(n):
                    if coalition >> player & 1:
                        continue
                    successor = coalition | (1 << player)
                    self.inequalities.append(
                        {
                            self.column(node, coalition): 1.0,
                            self.column(node, successor): -1.0,
                        }
                    )
                    self.inequality_rhs.append(0.0)
        for lower, upper, protected in topology.arcs:
            for coalition in range(self.coalition_count):
                row = {
                    self.column(lower, coalition): 1.0,
                    self.column(upper, coalition): -1.0,
                }
                if coalition & protected != protected:
                    self.equalities.append(row)
                    self.equality_rhs.append(0.0)
                else:
                    self.inequalities.append(row)
                    self.inequality_rhs.append(0.0)
        self.a_ub = sparse_matrix(self.inequalities, self.variable_count)
        self.b_ub = np.asarray(self.inequality_rhs)
        self.a_eq = sparse_matrix(self.equalities, self.variable_count)
        self.b_eq = np.asarray(self.equality_rhs)

    def column(self, node: int, coalition: int) -> int:
        return node * self.coalition_count + coalition

    def maximize(self, objective: np.ndarray) -> np.ndarray | None:
        result = linprog(
            -objective,
            A_ub=self.a_ub,
            b_ub=self.b_ub,
            A_eq=self.a_eq,
            b_eq=self.b_eq,
            bounds=[(None, None)] * self.variable_count,
            method="highs",
        )
        return result.x if result.success else None


def lower_expectation_dual(
    game: np.ndarray, divergence: tuple[int, ...], n: int
) -> tuple[float, np.ndarray, float] | None:
    grand = (1 << n) - 1
    proper = list(range(1, grand))
    beta_index = len(proper)
    objective = np.zeros(beta_index + 1)
    for index, coalition in enumerate(proper):
        objective[index] = -float(game[coalition])
    objective[beta_index] = 1.0
    equalities = np.zeros((n, beta_index + 1))
    for player in range(n):
        for index, coalition in enumerate(proper):
            if coalition >> player & 1:
                equalities[player, index] = 1.0
        equalities[player, beta_index] = -1.0
    result = linprog(
        objective,
        A_eq=equalities,
        b_eq=np.asarray(divergence, dtype=float),
        bounds=[(0.0, None)] * beta_index + [(None, None)],
        method="highs",
    )
    if not result.success:
        return None
    return -float(result.fun), result.x[:beta_index], float(
        result.x[beta_index]
    )


def random_topology(
    n: int, node_count: int, rng: random.Random
) -> Topology | None:
    arc_players: dict[tuple[int, int], int] = {}
    divergence = [[0] * n for _ in range(node_count)]
    for player in range(n):
        arc_count = rng.randint(1, max(2, node_count - 1))
        chosen = rng.sample(
            [
                (lower, upper)
                for lower in range(node_count)
                for upper in range(lower + 1, node_count)
            ],
            min(arc_count, node_count * (node_count - 1) // 2),
        )
        for lower, upper in chosen:
            weight = rng.randint(1, 3)
            arc_players[lower, upper] = (
                arc_players.get((lower, upper), 0) | 1 << player
            )
            divergence[lower][player] += weight
            divergence[upper][player] -= weight
    if any(not any(vector) for vector in divergence):
        return None
    if not any(len(set(vector)) >= 3 for vector in divergence):
        return None
    undirected = [set() for _ in range(node_count)]
    for lower, upper in arc_players:
        undirected[lower].add(upper)
        undirected[upper].add(lower)
    seen = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for neighbor in undirected[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    if len(seen) != node_count:
        return None
    return Topology(
        node_count,
        tuple(
            (lower, upper, protected)
            for (lower, upper), protected in sorted(arc_players.items())
        ),
        tuple(tuple(vector) for vector in divergence),
    )


def four_cycle_topology(weights: tuple[int, int, int, int], n: int) -> Topology:
    if n < 4:
        raise ValueError("the four-cycle topology needs at least four players")
    divergence = [[0] * n for _ in range(4)]
    arcs = ((0, 2, 1), (0, 3, 2), (1, 2, 4), (1, 3, 8))
    for (lower, upper, protected), weight in zip(arcs, weights, strict=True):
        player = protected.bit_length() - 1
        divergence[lower][player] += weight
        divergence[upper][player] -= weight
    return Topology(4, arcs, tuple(tuple(vector) for vector in divergence))


def optimize_topology(
    topology: Topology,
    n: int,
    rng: random.Random,
    starts: int,
    iterations: int,
) -> dict[str, Any] | None:
    cone = GameCone(n, topology)
    grand = (1 << n) - 1
    proper = list(range(1, grand))
    best: dict[str, Any] | None = None
    initial_points: list[np.ndarray] = []
    if n == 5:
        permutations = list(itertools.permutations(range(n)))
        rng.shuffle(permutations)
        for permutation in permutations[: min(12, starts)]:
            first = (3, 1, 1, 2, 2)
            second = (2, 3, 2, 1, 1)
            game = np.zeros(grand + 1)
            for coalition in range(grand + 1):
                game[coalition] = min(
                    sum(
                        first[permutation[player]]
                        for player in range(n)
                        if coalition >> player & 1
                    ),
                    sum(
                        second[permutation[player]]
                        for player in range(n)
                        if coalition >> player & 1
                    ),
                ) / 9.0
            initial_points.append(np.tile(game, topology.node_count))
    for _ in range(starts):
        random_objective = np.asarray(
            [rng.uniform(-1.0, 1.0) for _ in range(cone.variable_count)]
        )
        point = cone.maximize(random_objective)
        if point is not None:
            initial_points.append(point)
    for point in initial_points:
        values: list[float] = []
        for _iteration in range(iterations):
            objective = np.zeros(cone.variable_count)
            betas: list[float] = []
            values = []
            valid = True
            for node, divergence in enumerate(topology.divergence):
                game = point[
                    node * (grand + 1) : (node + 1) * (grand + 1)
                ]
                dual = lower_expectation_dual(game, divergence, n)
                if dual is None:
                    valid = False
                    break
                value, weights, beta = dual
                values.append(value)
                betas.append(beta)
                for index, coalition in enumerate(proper):
                    objective[cone.column(node, coalition)] += weights[index]
            if not valid:
                break
            updated = cone.maximize(objective)
            if updated is None:
                break
            point = updated
        final_values = []
        for node, divergence in enumerate(topology.divergence):
            game = point[
                node * (grand + 1) : (node + 1) * (grand + 1)
            ]
            dual = lower_expectation_dual(game, divergence, n)
            if dual is None:
                final_values = []
                break
            final_values.append(dual[0])
        if not final_values:
            continue
        total = sum(final_values)
        candidate = {
            "objective": total,
            "node_lower_expectations": final_values,
            "games": [
                [
                    float(point[node * (grand + 1) + coalition])
                    for coalition in range(grand + 1)
                ]
                for node in range(topology.node_count)
            ],
        }
        if best is None or total > float(best["objective"]):
            best = candidate
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--nodes", type=int, default=5)
    parser.add_argument("--samples", type=int, default=1_000)
    parser.add_argument("--starts", type=int, default=4)
    parser.add_argument("--iterations", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--canonical-four-cycle-max-weight", type=int)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    tested = 0
    positive = 0
    best: list[dict[str, Any]] = []
    if args.canonical_four_cycle_max_weight is not None:
        candidates: Any = (
            four_cycle_topology(weights, args.n)
            for weights in itertools.product(
                range(1, args.canonical_four_cycle_max_weight + 1),
                repeat=4,
            )
        )
        sample_target = args.canonical_four_cycle_max_weight**4
    else:
        candidates = (
            random_topology(args.n, args.nodes, rng)
            for _ in itertools.count()
        )
        sample_target = args.samples
    while tested < sample_target:
        topology = next(candidates)
        if topology is None:
            continue
        tested += 1
        optimum = optimize_topology(
            topology,
            args.n,
            rng,
            args.starts,
            args.iterations,
        )
        if optimum is None:
            continue
        row = {
            "arcs": [list(arc) for arc in topology.arcs],
            "divergence": [list(vector) for vector in topology.divergence],
            **optimum,
        }
        if float(row["objective"]) > 1e-8:
            positive += 1
        best.append(row)
        best.sort(key=lambda candidate: float(candidate["objective"]), reverse=True)
        del best[args.top :]
    payload = {
        "status": "mixed_terminal_coordinate_ascent_float",
        "n": args.n,
        "node_count": (
            4
            if args.canonical_four_cycle_max_weight is not None
            else args.nodes
        ),
        "tested_topologies": tested,
        "positive_topologies": positive,
        "best": best,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
