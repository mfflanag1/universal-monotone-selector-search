#!/usr/bin/env python3
"""Certify the five-player 2-1-2 nested lower-expectation identity."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import z3
from scipy.optimize import linprog
from scipy.optimize._highspy._core import (
    HighsBasisStatus,
    HighsLp,
    MatrixFormat,
    _Highs,
    kHighsInf,
)
from scipy.sparse import csc_matrix


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
sys.path.insert(0, str(PROJECT_SRC))

from n5_facet_search import load_facets

from n5_mixed_divergence_premium import choquet_representation
from n5_mixed_terminal_branch_search import extreme_representations


F = Fraction


def exact_z3_cone_solution(
    target: tuple[F, ...], generators: list[tuple[F, ...]]
) -> list[tuple[int, F]] | None:
    """Decide rational cone membership exactly with linear arithmetic."""
    weights = [z3.Real(f"cone_weight_{index}") for index in range(len(generators))]
    solver = z3.Solver()
    solver.set(timeout=300_000)
    solver.add(*(weight >= 0 for weight in weights))

    def rational(value: F):
        return z3.Q(value.numerator, value.denominator)

    for coordinate, expected in enumerate(target):
        terms = [
            weights[index] * rational(generator[coordinate])
            for index, generator in enumerate(generators)
            if generator[coordinate]
        ]
        solver.add((z3.Sum(terms) if terms else z3.RealVal(0)) == rational(expected))
    if solver.check() != z3.sat:
        return None
    model = solver.model()
    result: list[tuple[int, F]] = []
    for index, weight in enumerate(weights):
        value = model.eval(weight, model_completion=True)
        if not z3.is_rational_value(value):
            raise RuntimeError("non-rational value in linear-arithmetic model")
        exact = F(value.numerator_as_long(), value.denominator_as_long())
        if exact < 0:
            raise RuntimeError("negative coefficient in exact cone model")
        if exact:
            result.append((index, exact))
    reconstructed = [F(0)] * len(target)
    for index, weight in result:
        for coordinate, coefficient in enumerate(generators[index]):
            reconstructed[coordinate] += weight * coefficient
    if tuple(reconstructed) != target:
        raise RuntimeError("exact cone model failed reconstruction")
    return result


def exact_solution_on_support(
    target: tuple[F, ...],
    generators: list[tuple[F, ...]],
    support_indices: list[int],
    approximate_weights: list[float] | None = None,
    approximation_denominators: tuple[int, ...] | None = None,
) -> list[tuple[int, F]] | None:
    """Solve a small supported cone identity exactly by rational RREF."""
    if not support_indices:
        return [] if not any(target) else None
    column_count = len(support_indices)
    matrix = [
        [generators[index][row] for index in support_indices]
        + [target[row]]
        for row in range(len(target))
    ]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        selected = next(
            (
                row
                for row in range(pivot_row, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = (
            matrix[selected],
            matrix[pivot_row],
        )
        pivot = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [
                value - multiplier * pivot_value
                for value, pivot_value in zip(
                    matrix[row], matrix[pivot_row], strict=True
                )
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    if any(
        not any(row[:column_count]) and row[column_count]
        for row in matrix
    ):
        return None
    free_columns = [
        column
        for column in range(column_count)
        if column not in set(pivot_columns)
    ]
    denominators = (
        approximation_denominators
        if approximation_denominators is not None
        else (
            (1,)
            if approximate_weights is None
            else (1_000, 1_000_000, 1_000_000_000, 1_000_000_000_000)
        )
    )
    for denominator in denominators:
        weights = [F(0)] * column_count
        if approximate_weights is not None:
            for column in free_columns:
                weights[column] = F(
                    float(approximate_weights[column])
                ).limit_denominator(denominator)
        for row, column in enumerate(pivot_columns):
            weights[column] = matrix[row][column_count] - sum(
                (
                    matrix[row][free] * weights[free]
                    for free in free_columns
                ),
                F(0),
            )
        if any(weight < 0 for weight in weights):
            continue
        result = [
            (index, weight)
            for index, weight in zip(support_indices, weights, strict=True)
            if weight
        ]
        reconstructed = [F(0)] * len(target)
        for index, weight in result:
            for coordinate, coefficient in enumerate(generators[index]):
                reconstructed[coordinate] += weight * coefficient
        if tuple(reconstructed) == target:
            return result
    return None


def coefficient_vector(
    representation: tuple[F, tuple[tuple[int, F], ...]], n: int
) -> tuple[F, ...]:
    grand = (1 << n) - 1
    beta, terms = representation
    result = [F(0)] * (grand + 1)
    result[grand] = beta
    for coalition, weight in terms:
        result[coalition] += weight
    return tuple(result)


def exact_generators(n: int) -> tuple[list[tuple[F, ...]], list[tuple[Any, ...]]]:
    grand = (1 << n) - 1
    facets, _ = load_facets()
    generators = [tuple(F(value) for value in facet) for facet in facets]
    metadata: list[tuple[Any, ...]] = [
        ("exact_facet", index) for index in range(len(facets))
    ]
    for sign in (-1, 1):
        empty_row = [F(0)] * (grand + 1)
        empty_row[0] = F(sign)
        generators.append(tuple(empty_row))
        metadata.append(("empty_equality", sign))
    for coalition in range(grand + 1):
        for player in range(n):
            if coalition >> player & 1:
                continue
            successor = coalition | (1 << player)
            row = [F(0)] * (grand + 1)
            row[coalition] = F(-1)
            row[successor] = F(1)
            generators.append(tuple(row))
            metadata.append(
                ("game_monotonicity", coalition, successor)
            )
    return generators, metadata


def find_certificate(
    target: tuple[F, ...],
    generators: list[tuple[F, ...]],
    *,
    allow_z3: bool = True,
    alternate_basis_retries: int = 32,
) -> list[tuple[int, F]]:
    if not any(target):
        return []
    matrix = csc_matrix(
        np.asarray(
            [
                [float(generator[row]) for generator in generators]
                for row in range(len(target))
            ]
        )
    )
    objective = np.asarray(
        [1.0 + (index % 97) * 1e-7 for index in range(len(generators))]
    )
    solved = linprog(
        objective,
        A_eq=matrix,
        b_eq=np.asarray([float(value) for value in target]),
        bounds=[(0.0, None)] * len(generators),
        method="highs-ds",
        options={
            "dual_feasibility_tolerance": 1e-10,
            "primal_feasibility_tolerance": 1e-10,
        },
    )
    if not solved.success:
        raise RuntimeError("target is outside the generated dual cone")
    support = [
        (index, F(float(weight)).limit_denominator(1_000_000))
        for index, weight in enumerate(solved.x)
        if weight > 1e-9
    ]
    reconstructed = [F(0)] * len(target)
    for index, weight in support:
        for coordinate, coefficient in enumerate(generators[index]):
            reconstructed[coordinate] += weight * coefficient
    if tuple(reconstructed) != target:
        exact_support = None
        for tolerance in (1e-9, 1e-12, 1e-14, 0.0):
            support_indices = [
                index
                for index, weight in enumerate(solved.x)
                if weight > tolerance
            ]
            exact_support = exact_solution_on_support(
                target,
                generators,
                support_indices,
                [float(solved.x[index]) for index in support_indices],
            )
            if exact_support is not None:
                break
        if exact_support is None:
            lp = HighsLp()
            lp.num_col_ = len(generators)
            lp.num_row_ = len(target)
            lp.col_cost_ = objective
            lp.col_lower_ = np.zeros(len(generators))
            lp.col_upper_ = np.full(len(generators), kHighsInf)
            exact_rhs = np.asarray([float(value) for value in target])
            lp.row_lower_ = exact_rhs
            lp.row_upper_ = exact_rhs
            lp.a_matrix_.format_ = MatrixFormat.kColwise
            lp.a_matrix_.num_col_ = len(generators)
            lp.a_matrix_.num_row_ = len(target)
            lp.a_matrix_.start_ = matrix.indptr.astype(np.int32)
            lp.a_matrix_.index_ = matrix.indices.astype(np.int32)
            lp.a_matrix_.value_ = matrix.data
            highs = _Highs()
            highs.setOptionValue("output_flag", False)
            highs.setOptionValue("solver", "simplex")
            highs.setOptionValue("dual_feasibility_tolerance", 1e-10)
            highs.setOptionValue("primal_feasibility_tolerance", 1e-10)
            if highs.passModel(lp).name == "kOk" and highs.run().name == "kOk":
                basis = highs.getBasis()
                solution = highs.getSolution()
                support_indices = [
                    index
                    for index, status in enumerate(basis.col_status)
                    if status == HighsBasisStatus.kBasic
                    or solution.col_value[index] > 1e-12
                ]
                exact_support = exact_solution_on_support(
                    target,
                    generators,
                    support_indices,
                    [
                        float(solution.col_value[index])
                        for index in support_indices
                    ],
                )
        if exact_support is None:
            # Degenerate boundary points can yield a floating simplex basis
            # whose exact replay has a tiny negative basic coefficient.  Try
            # deterministic positive objective perturbations to expose an
            # alternate basis.  Any result is still accepted only after the
            # exact rational reconstruction in exact_solution_on_support.
            generator_indices = np.arange(len(generators), dtype=np.int64)
            for salt in range(1, alternate_basis_retries + 1):
                retry_objective = 1.0 + (
                    (
                        generator_indices * (2 * salt + 1)
                        + salt * salt
                    )
                    % 1009
                ) * 1e-6
                retried = linprog(
                    retry_objective,
                    A_eq=matrix,
                    b_eq=np.asarray([float(value) for value in target]),
                    bounds=[(0.0, None)] * len(generators),
                    method="highs-ds",
                    options={
                        "dual_feasibility_tolerance": 1e-10,
                        "primal_feasibility_tolerance": 1e-10,
                    },
                )
                if not retried.success:
                    continue
                for tolerance in (1e-9, 1e-12, 1e-14, 0.0):
                    support_indices = [
                        index
                        for index, weight in enumerate(retried.x)
                        if weight > tolerance
                    ]
                    exact_support = exact_solution_on_support(
                        target,
                        generators,
                        support_indices,
                        [
                            float(retried.x[index])
                            for index in support_indices
                        ],
                    )
                    if exact_support is not None:
                        break
                if exact_support is not None:
                    print(
                        json.dumps(
                            {
                                "exact_cone_alternate_basis_salt": salt,
                                "exact_cone_alternate_basis_support": len(
                                    exact_support
                                ),
                            }
                        ),
                        flush=True,
                    )
                    break
        if exact_support is None and allow_z3:
            exact_support = exact_z3_cone_solution(target, generators)
        if exact_support is None:
            raise RuntimeError(
                "floating support did not admit an exact nonnegative "
                "rational reconstruction"
            )
        support = exact_support
    if any(weight < 0 for _, weight in support):
        raise RuntimeError("certificate has a negative cone coefficient")
    return support


def generate(output: Path) -> None:
    n = 5
    divergence = (F(-1), F(-1), F(0), F(1), F(1))
    choquet = choquet_representation(divergence, n)
    choquet_vector = coefficient_vector(choquet, n)
    generators, metadata = exact_generators(n)
    certificates = []
    for branch, representation in enumerate(
        extreme_representations(list(divergence), n)
    ):
        branch_vector = coefficient_vector(representation, n)
        target = tuple(
            choquet_vector[index] - branch_vector[index]
            for index in range(1 << n)
        )
        support = find_certificate(target, generators)
        beta, terms = representation
        certificates.append(
            {
                "branch": branch,
                "beta": str(beta),
                "terms": [
                    [coalition, str(weight)] for coalition, weight in terms
                ],
                "dual_cone_support": [
                    {
                        "generator": list(metadata[index]),
                        "weight": str(weight),
                    }
                    for index, weight in support
                ],
            }
        )
    payload = {
        "status": "n5_middle_split_identity_exact_certificate",
        "n": n,
        "divergence": [str(value) for value in divergence],
        "choquet_beta": str(choquet[0]),
        "choquet_terms": [
            [coalition, str(weight)] for coalition, weight in choquet[1]
        ],
        "branch_count": len(certificates),
        "certificates": certificates,
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")


def verify(archive: Path) -> None:
    payload = json.loads(archive.read_text())
    if payload["status"] != "n5_middle_split_identity_exact_certificate":
        raise RuntimeError("unexpected archive status")
    n = int(payload["n"])
    divergence = tuple(F(value) for value in payload["divergence"])
    expected_representations = extreme_representations(list(divergence), n)
    if int(payload["branch_count"]) != len(expected_representations):
        raise RuntimeError("branch count mismatch")
    generators, metadata = exact_generators(n)
    generator_lookup = {
        tuple(row): index for index, row in enumerate(metadata)
    }
    choquet = (
        F(payload["choquet_beta"]),
        tuple(
            (int(coalition), F(weight))
            for coalition, weight in payload["choquet_terms"]
        ),
    )
    if choquet != choquet_representation(divergence, n):
        raise RuntimeError("Choquet representation mismatch")
    choquet_vector = coefficient_vector(choquet, n)
    seen = set()
    for entry in payload["certificates"]:
        branch = int(entry["branch"])
        if branch in seen or not 0 <= branch < len(expected_representations):
            raise RuntimeError("invalid or duplicate branch")
        seen.add(branch)
        representation = (
            F(entry["beta"]),
            tuple(
                (int(coalition), F(weight))
                for coalition, weight in entry["terms"]
            ),
        )
        if representation != expected_representations[branch]:
            raise RuntimeError("branch representation mismatch")
        branch_vector = coefficient_vector(representation, n)
        target = tuple(
            choquet_vector[index] - branch_vector[index]
            for index in range(1 << n)
        )
        reconstructed = [F(0)] * (1 << n)
        for active in entry["dual_cone_support"]:
            generator = tuple(active["generator"])
            if generator not in generator_lookup:
                raise RuntimeError("unknown dual-cone generator")
            weight = F(active["weight"])
            if weight < 0:
                raise RuntimeError("negative dual-cone weight")
            for coordinate, coefficient in enumerate(
                generators[generator_lookup[generator]]
            ):
                reconstructed[coordinate] += weight * coefficient
        if tuple(reconstructed) != target:
            raise RuntimeError("exact cone identity failed")
    if len(seen) != len(expected_representations):
        raise RuntimeError("archive omits a lower-expectation branch")
    print(
        "PASS: exact five-player 2-1-2 identity, "
        f"{len(seen)} lower-expectation branches"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--generate", action="store_true")
    args = parser.parse_args()
    if args.generate:
        generate(args.archive)
    verify(args.archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
