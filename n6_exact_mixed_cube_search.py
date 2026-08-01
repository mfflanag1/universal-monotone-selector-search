#!/usr/bin/env python3
"""Adversarial exact-cone search on six-player mixed grand/proper subcomplexes."""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

sys.path.append(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)

from family_search import family_constraints
from n6_mixed_terminal_search import load_facets


N = 6
GRAND = 63


def sparse(rows: list[dict[int, float]], columns: int):
    rr, cc, vv = [], [], []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                rr.append(row_index)
                cc.append(column)
                vv.append(value)
    return coo_matrix((vv, (rr, cc)), shape=(len(rows), columns)).tocsr()


class ExactTopology:
    def __init__(self, edges: list[tuple[int, int, int]], bumps: list[float], facets):
        self.edges = edges
        self.game_count = 1 + max(max(a, b) for a, b, _ in edges)
        self.variable_count = self.game_count * GRAND
        ub, ub_rhs, eq, eq_rhs = [], [], [], []
        eq.append({self.column(0, GRAND): 1.0})
        eq_rhs.append(1.0)
        for (lower, upper, changed), bump in zip(edges, bumps, strict=True):
            for coalition in range(1, GRAND + 1):
                eq.append({self.column(upper, coalition): 1.0, self.column(lower, coalition): -1.0})
                eq_rhs.append(bump if coalition == changed else 0.0)
        for node in range(self.game_count):
            for facet in facets:
                ub.append({self.column(node, coalition): -float(coefficient) for coalition, coefficient in enumerate(facet) if coalition and coefficient})
                ub_rhs.append(0.0)
            for coalition in range(GRAND):
                for player in range(N):
                    if coalition >> player & 1:
                        continue
                    successor = coalition | (1 << player)
                    row = {self.column(node, successor): -1.0}
                    if coalition:
                        row[self.column(node, coalition)] = 1.0
                    ub.append(row)
                    ub_rhs.append(0.0)
        self.a_ub = sparse(ub, self.variable_count)
        self.b_ub = np.asarray(ub_rhs)
        self.a_eq = sparse(eq, self.variable_count)
        self.b_eq = np.asarray(eq_rhs)

    def column(self, node: int, coalition: int) -> int:
        return node * GRAND + coalition - 1

    def solve(self, objective: np.ndarray):
        result = linprog(objective, A_ub=self.a_ub, b_ub=self.b_ub, A_eq=self.a_eq, b_eq=self.b_eq, bounds=[(None, None)] * self.variable_count, method="highs")
        if not result.success:
            return None
        return [[0.0] + [float(result.x[self.column(node, coalition)]) for coalition in range(1, GRAND + 1)] for node in range(self.game_count)]


def margin_dual(games, edges):
    inequalities, rhs, equalities, eq_rhs, metadata = family_constraints(games, edges, N)
    variables = len(games) * N + 1
    a_ub = []
    for row, tag in zip(inequalities, metadata, strict=True):
        extended = [float(value) for value in row] + [1.0 if tag[0] == "monotonicity" else 0.0]
        a_ub.append(extended)
    a_eq = [[float(value) for value in row] + [0.0] for row in equalities]
    objective = np.zeros(variables)
    objective[-1] = -1.0
    result = linprog(objective, A_ub=a_ub, b_ub=np.asarray(rhs, float), A_eq=a_eq, b_eq=np.asarray(eq_rhs, float), bounds=[(None, None)] * variables, method="highs")
    if not result.success:
        raise RuntimeError(result.message)
    return float(result.x[-1]), result.ineqlin.marginals, result.eqlin.marginals, metadata


def flow_objective(model, iw, ew, metadata):
    objective = np.zeros(model.variable_count)
    for weight, tag in zip(iw, metadata, strict=True):
        if weight and tag[0] == "core":
            objective[model.column(int(tag[1]), int(tag[2]))] += weight
    for node, weight in enumerate(ew):
        objective[model.column(node, GRAND)] -= weight
    return objective


def optimize(edges, bumps, facets, starts, iterations, seed):
    model = ExactTopology(edges, bumps, facets)
    rng = random.Random(seed)
    best_margin, best_games = float("inf"), None
    for _ in range(starts):
        objective = np.asarray([rng.uniform(-1, 1) for _ in range(model.variable_count)])
        games = model.solve(objective)
        if games is None:
            continue
        for _ in range(iterations):
            margin, iw, ew, metadata = margin_dual(games, edges)
            if margin < best_margin:
                best_margin, best_games = margin, games
            updated = model.solve(flow_objective(model, iw, ew, metadata))
            if updated is None:
                break
            new_margin = margin_dual(updated, edges)[0]
            games = updated
            if new_margin >= margin - 1e-9:
                if new_margin < best_margin:
                    best_margin, best_games = new_margin, games
                break
    return best_margin, best_games


def states_and_edges(dimension, node_count, directions, rng):
    states = {0}
    while len(states) < node_count:
        parent = rng.choice(tuple(states))
        choices = [bit for bit in range(dimension) if not parent >> bit & 1]
        if not choices:
            continue
        states.add(parent | (1 << rng.choice(choices)))
    ordered = sorted(states)
    index = {state: i for i, state in enumerate(ordered)}
    edges = []
    for state in ordered:
        for bit, coalition in enumerate(directions):
            child = state | (1 << bit)
            if child != state and child in states:
                edges.append((index[state], index[child], coalition))
    return edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=20)
    parser.add_argument("--dimension", type=int, default=7)
    parser.add_argument("--nodes", type=int, default=24)
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--bump", type=float, default=0.005)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument(
        "--directions",
        help="comma-separated fixed coalition masks (overrides --dimension)",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facets = load_facets()
    rng = random.Random(args.seed)
    rows, best = [], None
    for case in range(args.cases):
        if args.directions:
            directions = [int(value) for value in args.directions.split(",")]
            if not directions or any(value < 1 or value > GRAND for value in directions):
                raise ValueError("directions must be coalition masks between 1 and 63")
        else:
            directions = rng.sample(range(1, GRAND), args.dimension - 1) + [GRAND]
            rng.shuffle(directions)
        edges = states_and_edges(len(directions), args.nodes, directions, rng)
        if not args.directions and GRAND not in {coalition for _, _, coalition in edges}:
            continue
        by_direction = {coalition: args.bump * rng.randint(1, 6) / 3 for coalition in directions}
        bumps = [by_direction[coalition] for _, _, coalition in edges]
        margin, games = optimize(edges, bumps, facets, args.starts, args.iterations, args.seed + case)
        row = {"case": case, "directions": directions, "nodes": args.nodes, "edges": len(edges), "margin": margin}
        rows.append(row)
        if best is None or margin < best["margin"]:
            best = row | {"edge_list": edges, "bumps": bumps, "games": games}
            print(json.dumps({"best": row}), flush=True)
        if margin < -1e-8:
            break
    payload = {"status": "negative_found" if best and best["margin"] < -1e-8 else "no_negative_found", "rows": rows, "best": best}
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
