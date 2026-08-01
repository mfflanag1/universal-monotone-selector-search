#!/usr/bin/env python3
"""Mixed-divergence terminal search on the complete six-player exact cone."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import random
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog

from n5_mixed_terminal_search import (
    Topology,
    lower_expectation_dual,
    random_topology,
    sparse_matrix,
)


N = 6
GRAND = (1 << N) - 1
DEFAULT_CATALOGUE = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/data/n6_exact_facets.csv"
)


def label_mask(label: str) -> int:
    return sum(1 << (ord(character) - ord("a")) for character in label)


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(
        1 << permutation[player]
        for player in range(N)
        if mask >> player & 1
    )


def load_facets(path: Path = DEFAULT_CATALOGUE) -> list[tuple[int, ...]]:
    rows = list(csv.reader(path.open(), delimiter=";"))
    masks = [label_mask(label) if label else 0 for label in rows[0]]
    representatives = []
    for raw in rows[1:]:
        coefficients = [0] * (1 << N)
        for column, text in enumerate(raw):
            if column:
                coefficients[masks[column]] = int(text)
        representatives.append(tuple(coefficients))
    orbit = set()
    for representative in representatives:
        for permutation in itertools.permutations(range(N)):
            coefficients = [0] * (1 << N)
            for coalition, coefficient in enumerate(representative):
                if coefficient:
                    coefficients[permute_mask(coalition, permutation)] = coefficient
            orbit.add(tuple(coefficients))
    if len(orbit) != 7006:
        raise RuntimeError(f"expected 7006 facets, got {len(orbit)}")
    return list(orbit)


class GameCone6:
    def __init__(self, topology: Topology, facets: list[tuple[int, ...]]) -> None:
        self.node_count = topology.node_count
        self.coalition_count = 1 << N
        self.variable_count = self.node_count * self.coalition_count
        inequalities: list[dict[int, float]] = []
        inequality_rhs: list[float] = []
        equalities: list[dict[int, float]] = []
        equality_rhs: list[float] = []
        for node in range(self.node_count):
            equalities.append({self.column(node, 0): 1.0})
            equality_rhs.append(0.0)
            equalities.append({self.column(node, GRAND): 1.0})
            equality_rhs.append(1.0)
            for facet in facets:
                inequalities.append(
                    {
                        self.column(node, coalition): -float(coefficient)
                        for coalition, coefficient in enumerate(facet)
                        if coefficient
                    }
                )
                inequality_rhs.append(0.0)
            for coalition in range(self.coalition_count):
                for player in range(N):
                    if coalition >> player & 1:
                        continue
                    successor = coalition | (1 << player)
                    inequalities.append(
                        {
                            self.column(node, coalition): 1.0,
                            self.column(node, successor): -1.0,
                        }
                    )
                    inequality_rhs.append(0.0)
        for lower, upper, protected in topology.arcs:
            for coalition in range(self.coalition_count):
                row = {
                    self.column(lower, coalition): 1.0,
                    self.column(upper, coalition): -1.0,
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

    def column(self, node: int, coalition: int) -> int:
        return node * self.coalition_count + coalition

    def maximize(self, objective: np.ndarray) -> np.ndarray | None:
        result = linprog(
            -objective,
            A_ub=self.a_ub,
            b_ub=self.b_ub,
            A_eq=self.a_eq,
            b_eq=self.b_eq,
            bounds=[(0.0, 1.0)] * self.variable_count,
            method="highs",
        )
        return result.x if result.success else None


def optimize(
    topology: Topology,
    facets: list[tuple[int, ...]],
    rng: random.Random,
    starts: int,
    iterations: int,
) -> dict[str, Any] | None:
    cone = GameCone6(topology, facets)
    best = None
    seed_points = (
        (4.0, 2.0, 5.0, 12.0, 19.0, 9.0),
        (17.0, 25.0, 2.0, 1.0, 6.0, 0.0),
        (10.0, 11.0, 2.0, 9.0, 13.0, 6.0),
        (9.5, 10.25, 2.25, 9.25, 13.5, 6.25),
    )
    seed_game = np.asarray(
        [
            min(
                sum(point[player] for player in range(N) if coalition >> player & 1)
                for point in seed_points
            )
            / 51.0
            for coalition in range(1 << N)
        ]
    )
    initial_points = [np.tile(seed_game, topology.node_count)]
    for _ in range(starts):
        point = cone.maximize(
            np.asarray([rng.uniform(-1.0, 1.0) for _ in range(cone.variable_count)])
        )
        if point is None:
            continue
        initial_points.append(point)
    for point in initial_points:
        for _ in range(iterations):
            objective = np.zeros(cone.variable_count)
            for node, divergence in enumerate(topology.divergence):
                game = point[node * 64 : (node + 1) * 64]
                dual = lower_expectation_dual(game, divergence, N)
                if dual is None:
                    break
                _value, weights, _beta = dual
                for coalition, weight in enumerate(weights, start=1):
                    objective[cone.column(node, coalition)] += weight
            updated = cone.maximize(objective)
            if updated is None:
                break
            point = updated
        values = []
        for node, divergence in enumerate(topology.divergence):
            game = point[node * 64 : (node + 1) * 64]
            dual = lower_expectation_dual(game, divergence, N)
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
                "games": [
                    [float(value) for value in point[node * 64 : (node + 1) * 64]]
                    for node in range(topology.node_count)
                ],
            }
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=20)
    parser.add_argument("--nodes", type=int, default=4)
    parser.add_argument("--starts", type=int, default=4)
    parser.add_argument("--iterations", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facets = load_facets()
    rng = random.Random(args.seed)
    best = []
    tested = 0
    while tested < args.samples:
        topology = random_topology(N, args.nodes, rng)
        if topology is None:
            continue
        tested += 1
        optimum = optimize(topology, facets, rng, args.starts, args.iterations)
        if optimum is None:
            continue
        row = {
            "arcs": [list(arc) for arc in topology.arcs],
            "divergence": [list(vector) for vector in topology.divergence],
            **optimum,
        }
        best.append(row)
        best.sort(key=lambda item: float(item["objective"]), reverse=True)
        del best[10:]
        print(json.dumps({"tested": tested, "objective": row["objective"]}), flush=True)
        if float(row["objective"]) > 1e-8:
            break
    payload = {
        "status": (
            "positive_terminal_obstruction_found"
            if best and float(best[0]["objective"]) > 1e-8
            else "no_positive_terminal_obstruction_found"
        ),
        "facet_count": len(facets),
        "tested": tested,
        "configuration": vars(args) | {"output": str(args.output)},
        "best": best,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
