#!/usr/bin/env python3
"""Exhaust the two-game cross-cosingleton lower-expectation branches.

The Johnson-square objective is bounded above by this two-terminal swap
problem.  Each lower expectation is the maximum of 397 affine core-dual
branches.  This script fixes a pair of branches, warm-starts the same LP,
and records the largest branch-pair optimum.  It is a discovery pass; a
nonpositive floating result still needs exact dual reconstruction.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize._highspy._core import (
    HighsLp,
    MatrixFormat,
    ObjSense,
    _Highs,
    kHighsInf,
)
from scipy.sparse import vstack

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_mixed_terminal_branch_search import extreme_representations
from n5_mixed_terminal_search import (
    GameCone,
    cross_cosingleton_swap_topology,
)


F = Fraction


def branch_vector(
    representation: tuple[F, tuple[tuple[int, F], ...]],
    coalition_count: int,
) -> np.ndarray:
    beta, terms = representation
    vector = np.zeros(coalition_count)
    vector[-1] = float(beta)
    for coalition, weight in terms:
        vector[coalition] += float(weight)
    return vector


def exact_branch_vector(
    representation: tuple[F, tuple[tuple[int, F], ...]],
    coalition_count: int,
) -> list[F]:
    beta, terms = representation
    vector = [F(0)] * coalition_count
    vector[-1] = beta
    for coalition, weight in terms:
        vector[coalition] += weight
    return vector


def audit_exact_dual(
    row_dual: np.ndarray,
    exact_rows: list[dict[int, F]],
    exact_rhs: list[F],
    inequality_count: int,
    exact_cost: list[F],
) -> int | None:
    support = [
        index for index, value in enumerate(row_dual) if abs(value) > 1e-9
    ]
    for denominator in (1_000, 10_000, 100_000, 1_000_000, 100_000_000):
        weights = {
            index: F(float(row_dual[index])).limit_denominator(denominator)
            for index in support
        }
        if any(
            weight > 0
            for index, weight in weights.items()
            if index < inequality_count
        ):
            continue
        stationarity = [F(0)] * len(exact_cost)
        dual_objective = F(0)
        for index, weight in weights.items():
            for column, coefficient in exact_rows[index].items():
                stationarity[column] += weight * coefficient
            dual_objective += weight * exact_rhs[index]
        if stationarity == exact_cost and dual_objective == 0:
            return denominator
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int)
    parser.add_argument("--progress-every", type=int, default=10_000)
    parser.add_argument("--exact-audit", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    topology = cross_cosingleton_swap_topology(5)
    cone = GameCone(5, topology)
    representations = [
        extreme_representations([F(value) for value in divergence], 5)
        for divergence in topology.divergence
    ]
    vectors = [
        [branch_vector(rep, cone.coalition_count) for rep in choices]
        for choices in representations
    ]
    exact_vectors = [
        [exact_branch_vector(rep, cone.coalition_count) for rep in choices]
        for choices in representations
    ]

    matrix = vstack((cone.a_ub, cone.a_eq)).tocsc()
    exact_rows = [
        {column: F(value).limit_denominator() for column, value in row.items()}
        for row in cone.inequalities + cone.equalities
    ]
    exact_rhs = [
        F(value).limit_denominator()
        for value in cone.inequality_rhs + cone.equality_rhs
    ]
    lp = HighsLp()
    lp.num_col_ = cone.variable_count
    lp.num_row_ = matrix.shape[0]
    lp.col_cost_ = np.zeros(cone.variable_count)
    lp.col_lower_ = np.full(cone.variable_count, -kHighsInf)
    lp.col_upper_ = np.full(cone.variable_count, kHighsInf)
    lp.row_lower_ = np.concatenate(
        (np.full(cone.a_ub.shape[0], -kHighsInf), cone.b_eq)
    )
    lp.row_upper_ = np.concatenate((cone.b_ub, cone.b_eq))
    lp.a_matrix_.format_ = MatrixFormat.kColwise
    lp.a_matrix_.num_col_ = cone.variable_count
    lp.a_matrix_.num_row_ = matrix.shape[0]
    lp.a_matrix_.start_ = matrix.indptr.astype(np.int32)
    lp.a_matrix_.index_ = matrix.indices.astype(np.int32)
    lp.a_matrix_.value_ = matrix.data
    lp.sense_ = ObjSense.kMinimize

    highs = _Highs()
    highs.setOptionValue("output_flag", False)
    highs.setOptionValue("solver", "simplex")
    highs.setOptionValue("simplex_strategy", 1)
    if highs.passModel(lp).name != "kOk" or highs.run().name != "kOk":
        raise RuntimeError("HiGHS failed to initialize the branch sweep")

    indices = np.arange(cone.variable_count, dtype=np.int32)
    maximum = -float("inf")
    worst_pair: tuple[int, int] | None = None
    positive_count = 0
    exact_certificate_count = 0
    exact_reconstruction_max_denominator = 1
    exact_failures: list[list[int]] = []
    solved = 0
    started = time.monotonic()
    total = len(vectors[0]) * len(vectors[1])
    if args.limit is not None:
        total = min(total, args.limit)
    for left, right in itertools.product(
        range(len(vectors[0])), range(len(vectors[1]))
    ):
        if solved >= total:
            break
        objective = -np.concatenate((vectors[0][left], vectors[1][right]))
        if highs.changeColsCost(
            cone.variable_count, indices, objective
        ).name != "kOk":
            raise RuntimeError("HiGHS rejected an objective update")
        if highs.run().name != "kOk":
            raise RuntimeError(f"HiGHS failed at branch pair {(left, right)}")
        value = -float(highs.getObjectiveValue())
        solved += 1
        if args.exact_audit:
            exact_cost = [
                -value
                for value in exact_vectors[0][left] + exact_vectors[1][right]
            ]
            solution = highs.getSolution()
            denominator = audit_exact_dual(
                np.asarray(solution.row_dual),
                exact_rows,
                exact_rhs,
                cone.a_ub.shape[0],
                exact_cost,
            )
            if denominator is None:
                exact_failures.append([left, right])
            else:
                exact_certificate_count += 1
                exact_reconstruction_max_denominator = max(
                    exact_reconstruction_max_denominator, denominator
                )
        if value > maximum:
            maximum = value
            worst_pair = (left, right)
        if value > 1e-8:
            positive_count += 1
        if args.progress_every and solved % args.progress_every == 0:
            elapsed = time.monotonic() - started
            print(
                json.dumps(
                    {
                        "solved": solved,
                        "total": total,
                        "maximum": maximum,
                        "worst_pair": worst_pair,
                        "positive_count": positive_count,
                        "exact_certificate_count": exact_certificate_count,
                        "exact_failure_count": len(exact_failures),
                        "elapsed_seconds": elapsed,
                    }
                ),
                flush=True,
            )

    payload = {
        "status": (
            "cross_cosingleton_branch_sweep_exact"
            if args.exact_audit
            and solved == len(vectors[0]) * len(vectors[1])
            and exact_certificate_count == solved
            else "cross_cosingleton_branch_sweep_float"
        ),
        "representation_counts": [len(choices) for choices in representations],
        "solved_branch_pairs": solved,
        "total_branch_pairs": len(vectors[0]) * len(vectors[1]),
        "maximum_float": maximum,
        "worst_pair": list(worst_pair) if worst_pair is not None else None,
        "positive_count_at_1e-8": positive_count,
        "exact_audit_requested": args.exact_audit,
        "exact_certificate_count": exact_certificate_count,
        "exact_failure_count": len(exact_failures),
        "first_exact_failures": exact_failures[:20],
        "exact_reconstruction_max_denominator": (
            exact_reconstruction_max_denominator
        ),
        "exact_global_upper_bound": (
            "0"
            if args.exact_audit
            and solved == len(vectors[0]) * len(vectors[1])
            and exact_certificate_count == solved
            else None
        ),
        "elapsed_seconds": time.monotonic() - started,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if args.exact_audit and exact_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
