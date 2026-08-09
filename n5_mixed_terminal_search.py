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
    def __init__(
        self, n: int, topology: Topology, variable_grand: bool = False
    ) -> None:
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
            if variable_grand:
                self.inequalities.append(
                    {self.column(node, self.grand): 1.0}
                )
                self.inequality_rhs.append(1.0)
            else:
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


def johnson_square_topology(n: int) -> Topology:
    """The J(5,3) square making all four terminal premiums sharp."""
    if n != 5:
        raise ValueError("the Johnson-square topology is five-player")
    return block_square_topology((0b00111, 0b01011, 0b10101, 0b11001), n)


def block_square_topology(
    masks: tuple[int, int, int, int], n: int
) -> Topology:
    """A K2,2 protected-block cycle in row-major edge order."""
    grand = (1 << n) - 1
    if any(mask <= 0 or mask >= grand for mask in masks):
        raise ValueError("block-square masks must be nonempty and proper")
    arcs = tuple(
        (lower, upper, protected)
        for (lower, upper), protected in zip(
            ((0, 2), (0, 3), (1, 2), (1, 3)), masks, strict=True
        )
    )
    divergence = [[0] * n for _ in range(4)]
    for lower, upper, protected in arcs:
        for player in range(n):
            if protected >> player & 1:
                divergence[lower][player] += 1
                divergence[upper][player] -= 1
    return Topology(4, arcs, tuple(tuple(vector) for vector in divergence))


def bipartite_six_cycle_topology(
    masks: tuple[int, int, int, int, int, int], n: int
) -> Topology:
    """Alternating three-source/three-sink cycle with no universal hub."""
    return weighted_bipartite_six_cycle_topology(
        masks, (1, 1, 1, 1, 1, 1), n
    )


def weighted_bipartite_six_cycle_topology(
    masks: tuple[int, int, int, int, int, int],
    weights: tuple[int, int, int, int, int, int],
    n: int,
) -> Topology:
    """Weighted alternating three-source/three-sink protected cycle."""
    grand = (1 << n) - 1
    if any(mask <= 0 or mask >= grand for mask in masks):
        raise ValueError("six-cycle masks must be nonempty and proper")
    if any(weight <= 0 for weight in weights):
        raise ValueError("six-cycle weights must be positive")
    endpoints = ((0, 3), (1, 3), (1, 4), (2, 4), (2, 5), (0, 5))
    arcs = tuple(
        (lower, upper, protected)
        for (lower, upper), protected in zip(endpoints, masks, strict=True)
    )
    divergence = [[0] * n for _ in range(6)]
    for (lower, upper, protected), weight in zip(arcs, weights, strict=True):
        for player in range(n):
            if protected >> player & 1:
                divergence[lower][player] += weight
                divergence[upper][player] -= weight
    return Topology(6, arcs, tuple(tuple(vector) for vector in divergence))


def complete_bipartite_topology(
    masks: tuple[int, ...],
    weights: tuple[int, ...],
    side: int,
    n: int,
) -> Topology:
    """Complete source/sink block flow in row-major source/sink order."""
    return complete_bipartite_rectangular_topology(
        masks, weights, side, side, n
    )


def complete_bipartite_rectangular_topology(
    masks: tuple[int, ...],
    weights: tuple[int, ...],
    source_count: int,
    sink_count: int,
    n: int,
) -> Topology:
    """Rectangular complete source/sink block flow in row-major order."""
    grand = (1 << n) - 1
    if len(masks) != source_count * sink_count or len(weights) != len(masks):
        raise ValueError("complete-bipartite data has the wrong entry count")
    if any(mask <= 0 or mask >= grand for mask in masks):
        raise ValueError("complete-bipartite masks must be nonempty and proper")
    if any(weight <= 0 for weight in weights):
        raise ValueError("complete-bipartite weights must be positive")
    arcs = tuple(
        (
            source,
            source_count + sink,
            masks[source * sink_count + sink],
        )
        for source in range(source_count)
        for sink in range(sink_count)
    )
    divergence = [[0] * n for _ in range(source_count + sink_count)]
    for (lower, upper, protected), weight in zip(arcs, weights, strict=True):
        for player in range(n):
            if protected >> player & 1:
                divergence[lower][player] += weight
                divergence[upper][player] -= weight
    return Topology(
        source_count + sink_count,
        arcs,
        tuple(tuple(vector) for vector in divergence),
    )


def cross_cosingleton_swap_topology(n: int) -> Topology:
    """Two-game interval-overlap reduction of the Johnson square.

    Node 0 may exceed node 1 only at N\\{3}; node 1 may exceed node 0
    only at N\\{2}.  The opposite divergences measure the two copies of
    coalitions 013 and 034.  A nonpositive optimum for this relaxation is
    sufficient to close the four-terminal Johnson square.
    """
    if n != 5:
        raise ValueError("the cross-cosingleton swap is five-player")
    grand = (1 << n) - 1
    omit_2 = grand ^ (1 << 2)
    omit_3 = grand ^ (1 << 3)
    positive = (2, 1, 0, 2, 1)
    return Topology(
        2,
        ((0, 1, omit_2), (1, 0, omit_3)),
        (positive, tuple(-value for value in positive)),
    )


def optimize_topology(
    topology: Topology,
    n: int,
    rng: random.Random,
    starts: int,
    iterations: int,
    variable_grand: bool = False,
) -> dict[str, Any] | None:
    cone = GameCone(n, topology, variable_grand=variable_grand)
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
    parser.add_argument("--johnson-square", action="store_true")
    parser.add_argument("--cross-cosingleton-swap", action="store_true")
    parser.add_argument("--block-square-masks", type=int, nargs=4)
    parser.add_argument("--six-cycle-masks", type=int, nargs=6)
    parser.add_argument("--six-cycle-weights", type=int, nargs=6)
    parser.add_argument("--k33-masks", type=int, nargs=9)
    parser.add_argument("--k33-weights", type=int, nargs=9)
    parser.add_argument("--k23-masks", type=int, nargs=6)
    parser.add_argument("--k23-weights", type=int, nargs=6)
    parser.add_argument("--variable-grand", action="store_true")
    args = parser.parse_args()
    if args.six_cycle_weights is not None and args.six_cycle_masks is None:
        parser.error("--six-cycle-weights requires --six-cycle-masks")
    if args.k33_weights is not None and args.k33_masks is None:
        parser.error("--k33-weights requires --k33-masks")
    if args.k23_weights is not None and args.k23_masks is None:
        parser.error("--k23-weights requires --k23-masks")
    rng = random.Random(args.seed)
    tested = 0
    positive = 0
    best: list[dict[str, Any]] = []
    if args.block_square_masks is not None:
        candidates = iter(
            (block_square_topology(tuple(args.block_square_masks), args.n),)
        )
        sample_target = 1
    elif args.six_cycle_masks is not None:
        weights = (
            tuple(args.six_cycle_weights)
            if args.six_cycle_weights is not None
            else (1, 1, 1, 1, 1, 1)
        )
        candidates = iter(
            (
                weighted_bipartite_six_cycle_topology(
                    tuple(args.six_cycle_masks), weights, args.n
                ),
            )
        )
        sample_target = 1
    elif args.k33_masks is not None:
        weights = (
            tuple(args.k33_weights)
            if args.k33_weights is not None
            else (1,) * 9
        )
        candidates = iter(
            (
                complete_bipartite_topology(
                    tuple(args.k33_masks), weights, 3, args.n
                ),
            )
        )
        sample_target = 1
    elif args.k23_masks is not None:
        weights = (
            tuple(args.k23_weights)
            if args.k23_weights is not None
            else (1,) * 6
        )
        candidates = iter(
            (
                complete_bipartite_rectangular_topology(
                    tuple(args.k23_masks), weights, 2, 3, args.n
                ),
            )
        )
        sample_target = 1
    elif args.cross_cosingleton_swap:
        candidates = iter((cross_cosingleton_swap_topology(args.n),))
        sample_target = 1
    elif args.johnson_square:
        candidates = iter((johnson_square_topology(args.n),))
        sample_target = 1
    elif args.canonical_four_cycle_max_weight is not None:
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
            args.variable_grand,
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
            2
            if args.cross_cosingleton_swap
            else 5
            if args.k23_masks is not None
            else 6
            if args.six_cycle_masks is not None or args.k33_masks is not None
            else 4
            if args.johnson_square
            or args.block_square_masks is not None
            or args.canonical_four_cycle_max_weight is not None
            else args.nodes
        ),
        "tested_topologies": tested,
        "positive_topologies": positive,
        "variable_grand": args.variable_grand,
        "best": best,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
