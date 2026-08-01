#!/usr/bin/env python3
"""Alternating screen for adverse flows on legal exact bowtie complexes."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from n5_actual_bowtie_milp import GRAND, LegalExactCone, bowtie
from n5_mixed_terminal_search import lower_expectation_dual


def matrices(cone: LegalExactCone):
    ub, ub_rhs, eq, eq_rhs = [], [], [], []
    for row, lower, upper in zip(cone.rows, cone.lower, cone.upper, strict=True):
        dense = np.zeros(cone.game_variables)
        for column, value in row.items():
            dense[column] = value
        if lower == upper:
            eq.append(dense)
            eq_rhs.append(lower)
        else:
            if np.isfinite(upper):
                ub.append(dense)
                ub_rhs.append(upper)
            if np.isfinite(lower):
                ub.append(-dense)
                ub_rhs.append(-lower)
    return tuple(np.asarray(value) for value in (ub, ub_rhs, eq, eq_rhs))


def solve_linear(cone, matrices_value, objective):
    ub, ub_rhs, eq, eq_rhs = matrices_value
    result = linprog(
        -objective,
        A_ub=ub,
        b_ub=ub_rhs,
        A_eq=eq,
        b_eq=eq_rhs,
        bounds=[(0.0, 1.0)] * cone.game_variables,
        method="highs",
    )
    return result.x if result.success else None


def flow_value(cone, topology, point):
    objective = np.zeros(cone.game_variables)
    values = []
    for node, divergence in enumerate(topology.divergence):
        game = point[node * (GRAND + 1) : (node + 1) * (GRAND + 1)]
        dual = lower_expectation_dual(game, divergence, 5)
        if dual is None:
            return None
        value, weights, _beta = dual
        values.append(value)
        for coalition, weight in enumerate(weights, start=1):
            objective[cone.column(node, coalition)] += weight
    return sum(values), values, objective


def optimize(topology, bump, starts, iterations, rng):
    cone = LegalExactCone(topology, bump)
    matrix_data = matrices(cone)
    best = None
    for _start in range(starts):
        objective = np.asarray(
            [rng.uniform(-1.0, 1.0) for _ in range(cone.game_variables)]
        )
        point = solve_linear(cone, matrix_data, objective)
        if point is None:
            continue
        for _iteration in range(iterations):
            evaluated = flow_value(cone, topology, point)
            if evaluated is None:
                break
            value, values, objective = evaluated
            if best is None or value > best["objective"]:
                best = {
                    "objective": value,
                    "node_lower_expectations": values,
                    "games": [
                        [
                            float(point[cone.column(node, coalition)])
                            for coalition in range(GRAND + 1)
                        ]
                        for node in range(topology.node_count)
                    ],
                }
            updated = solve_linear(cone, matrix_data, objective)
            if updated is None:
                break
            next_value = flow_value(cone, topology, updated)[0]
            point = updated
            if next_value <= value + 1e-9:
                break
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=100)
    parser.add_argument("--arms", type=int, default=2)
    parser.add_argument("--starts", type=int, default=5)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--minimum-bump", type=float, default=0.001)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pool = [mask for mask in range(1, GRAND) if 2 <= mask.bit_count() <= 4]
    rows = []
    best = None
    for case in range(args.cases):
        directions = tuple(tuple(rng.sample(pool, 2)) for _ in range(args.arms))
        topology = bowtie(directions, rng)
        result = optimize(
            topology, args.minimum_bump, args.starts, args.iterations, rng
        )
        if result is None:
            continue
        row = {
            "case": case,
            "directions": [list(pair) for pair in directions],
            "edges": [list(edge) for edge in topology.edges],
            "flows": [list(flow) for flow in topology.flows],
            "divergence": [list(value) for value in topology.divergence],
            **result,
        }
        rows.append({key: value for key, value in row.items() if key != "games"})
        if best is None or row["objective"] > best["objective"]:
            best = row
            print(
                json.dumps(
                    {
                        "case": case,
                        "objective": row["objective"],
                        "directions": row["directions"],
                    }
                ),
                flush=True,
            )
        if row["objective"] > 1e-8:
            break
    payload = {
        "status": (
            "positive_flow_found"
            if best is not None and best["objective"] > 1e-8
            else "none_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "best": best,
        "rows": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
