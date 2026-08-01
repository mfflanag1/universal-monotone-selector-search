#!/usr/bin/env python3
"""Globally optimize a protected-support upper-expectation drop on n=5."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_facet_search import load_facets

from n5_mixed_terminal_milp import lower_dual_vertices


N = 5
GRAND = (1 << N) - 1
GAME_VARIABLES = GRAND + 1
PAYOFF_OFFSET = GAME_VARIABLES
VARIABLE_COUNT = GAME_VARIABLES + N


def sparse(rows: list[dict[int, float]]):
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
        shape=(len(rows), VARIABLE_COUNT),
    ).tocsr()


def solve(support: int, weights: tuple[int, ...]) -> dict[str, object]:
    facets, _ = load_facets()
    objective_weights = tuple(
        weights[player] if support >> player & 1 else 0
        for player in range(N)
    )
    candidates = lower_dual_vertices(
        tuple(-weight for weight in objective_weights), N
    )
    inequalities: list[dict[int, float]] = []
    inequality_rhs: list[float] = []
    inequality_metadata: list[tuple[object, ...]] = []
    equalities: list[dict[int, float]] = [
        {0: 1.0},
        {GRAND: 1.0},
        {PAYOFF_OFFSET + player: 1.0 for player in range(N)},
    ]
    equality_rhs = [0.0, 1.0, 1.0]
    for facet_index, facet in enumerate(facets):
        inequalities.append(
            {
                coalition: -float(coefficient)
                for coalition, coefficient in enumerate(facet)
                if coefficient
            }
        )
        inequality_rhs.append(0.0)
        inequality_metadata.append(("exact_facet", facet_index))
    for coalition in range(GRAND + 1):
        for player in range(N):
            if coalition >> player & 1:
                continue
            inequalities.append(
                {coalition: 1.0, coalition | (1 << player): -1.0}
            )
            inequality_rhs.append(0.0)
            inequality_metadata.append(
                ("game_monotonicity", coalition, player)
            )
    for coalition in range(1, GRAND):
        if coalition & support == support:
            continue
        row = {coalition: 1.0}
        for player in range(N):
            if coalition >> player & 1:
                row[PAYOFF_OFFSET + player] = -1.0
        inequalities.append(row)
        inequality_rhs.append(0.0)
        inequality_metadata.append(("relaxed_core", coalition))
    a_ub = sparse(inequalities)
    b_ub = np.asarray(inequality_rhs)
    a_eq = sparse(equalities)
    b_eq = np.asarray(equality_rhs)
    bounds = [(0.0, 1.0)] * GAME_VARIABLES + [(None, None)] * N
    best: dict[str, object] | None = None
    completed = 0
    unbounded = 0
    for candidate_index, candidate in enumerate(candidates):
        # The optimized gap is b*y + E_v(-b).  The latter is the maximum
        # over these affine dual pieces, so the two maxima commute.
        objective = np.zeros(VARIABLE_COUNT)
        for player, weight in enumerate(objective_weights):
            objective[PAYOFF_OFFSET + player] = -float(weight)
        for coalition, coefficient in candidate["coefficients"].items():
            objective[coalition] = -float(coefficient)
        result = linprog(
            objective,
            A_ub=a_ub,
            b_ub=b_ub,
            A_eq=a_eq,
            b_eq=b_eq,
            bounds=bounds,
            method="highs",
        )
        if result.status == 3:
            unbounded += 1
            continue
        if not result.success:
            raise RuntimeError(
                f"candidate {candidate_index}: {result.message}"
            )
        completed += 1
        value = -float(result.fun) + float(candidate["constant"])
        if best is None or value > float(best["gap"]):
            active_dual = []
            for row_index, marginal in enumerate(result.ineqlin.marginals):
                multiplier = -float(marginal)
                if multiplier > 1e-9:
                    active_dual.append(
                        {
                            "row": row_index,
                            "metadata": list(inequality_metadata[row_index]),
                            "multiplier": multiplier,
                        }
                    )
            best = {
                "candidate_index": candidate_index,
                "gap": value,
                "candidate_constant": str(candidate["constant"]),
                "candidate_coefficients": {
                    str(coalition): str(coefficient)
                    for coalition, coefficient in candidate[
                        "coefficients"
                    ].items()
                },
                "game": [float(value) for value in result.x[:GAME_VARIABLES]],
                "relaxed_optimizer": [
                    float(value) for value in result.x[PAYOFF_OFFSET:]
                ],
                "active_inequality_dual": active_dual,
                "equality_dual": [
                    -float(value) for value in result.eqlin.marginals
                ],
            }
    return {
        "status": (
            "positive_global_gap_found"
            if best is not None and float(best["gap"]) > 1e-8
            else "global_gap_zero"
        ),
        "support": support,
        "weights": list(objective_weights),
        "candidate_count": len(candidates),
        "completed_candidates": completed,
        "unbounded_candidates": unbounded,
        "best": best,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--support", type=int, required=True)
    parser.add_argument("--weights", type=int, nargs=N, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.support <= 0 or args.support >= GRAND:
        raise ValueError("support must be nonempty and proper")
    if any(
        (args.weights[player] > 0) != bool(args.support >> player & 1)
        for player in range(N)
    ):
        raise ValueError("weights must be positive exactly on the support")
    payload = solve(args.support, tuple(args.weights))
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if payload["status"] == "positive_global_gap_found" else 0


if __name__ == "__main__":
    raise SystemExit(main())
