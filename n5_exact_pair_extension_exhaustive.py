#!/usr/bin/env python3
"""Exhaustively test protected core extension for five-player exact bumps.

For every bumped coalition and every extreme balanced collection, maximize
the Bondareva--Shapley violation of the upper core augmented by lower-point
bounds on protected players.  A positive optimum is exactly a lower core
point with no protected extension to the upper core.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


F = Fraction
N = 5
GRAND = (1 << N) - 1


def balanced_vertices() -> list[tuple[tuple[int, ...], tuple[F, ...]]]:
    incidence = {
        coalition: tuple((coalition >> player) & 1 for player in range(N))
        for coalition in range(1, GRAND + 1)
    }
    vertices: list[tuple[tuple[int, ...], tuple[F, ...]]] = []
    for width in range(1, N + 1):
        for support in itertools.combinations(range(1, GRAND + 1), width):
            matrix = np.asarray([incidence[c] for c in support], dtype=float).T
            if np.linalg.matrix_rank(matrix) != width:
                continue
            weights_float, *_ = np.linalg.lstsq(
                matrix, np.ones(N), rcond=None
            )
            if (
                np.min(weights_float) <= 1e-10
                or np.max(np.abs(matrix @ weights_float - 1.0)) > 1e-9
            ):
                continue
            weights = tuple(F(float(value)).limit_denominator(120) for value in weights_float)
            if all(
                sum(
                    weights[index] * incidence[coalition][player]
                    for index, coalition in enumerate(support)
                )
                == 1
                for player in range(N)
            ):
                vertices.append((support, weights))
    return vertices


def load_facets(project_src: Path) -> list[tuple[int, ...]]:
    sys.path.insert(0, str(project_src))
    from n5_facet_search import load_facets as project_load_facets

    facets, _ = project_load_facets()
    return facets


class PairModel:
    """LP for an exact pair and one arbitrary point in the lower core."""

    def __init__(self, coalition: int, facets: list[tuple[int, ...]]) -> None:
        self.coalition = coalition
        self.delta_index = GRAND
        self.x_start = GRAND + 1
        self.variable_count = GRAND + 1 + N
        rows: list[dict[int, float]] = []
        rhs: list[float] = []

        def v_index(mask: int) -> int:
            return mask - 1

        def add_ub(row: dict[int, float], bound: float) -> None:
            rows.append({key: value for key, value in row.items() if value})
            rhs.append(bound)

        # Exact-cone inequalities for v and w=v+delta*1_S.
        for facet in facets:
            add_ub(
                {
                    v_index(mask): -float(coefficient)
                    for mask, coefficient in enumerate(facet)
                    if mask and coefficient
                },
                0.0,
            )
            upper_row = {
                v_index(mask): -float(coefficient)
                for mask, coefficient in enumerate(facet)
                if mask and coefficient
            }
            if facet[coalition]:
                upper_row[self.delta_index] = -float(facet[coalition])
            add_ub(upper_row, 0.0)

        # Monotonicity for v and w.  This is the project's normalized exact
        # cone section and bounds every worth between zero and the grand worth.
        for lower in range(GRAND):
            for player in range(N):
                if lower >> player & 1:
                    continue
                upper = lower | (1 << player)
                lower_row: dict[int, float] = {}
                if lower:
                    lower_row[v_index(lower)] = 1.0
                lower_row[v_index(upper)] = -1.0
                add_ub(lower_row, 0.0)

                upper_row = dict(lower_row)
                if lower == coalition:
                    upper_row[self.delta_index] = (
                        upper_row.get(self.delta_index, 0.0) + 1.0
                    )
                if upper == coalition:
                    upper_row[self.delta_index] = (
                        upper_row.get(self.delta_index, 0.0) - 1.0
                    )
                add_ub(upper_row, 0.0)

        # x belongs to C(v).
        for mask in range(1, GRAND):
            row = {v_index(mask): 1.0}
            for player in range(N):
                if mask >> player & 1:
                    row[self.x_start + player] = -1.0
            add_ub(row, 0.0)

        row_indices: list[int] = []
        columns: list[int] = []
        values: list[float] = []
        for row_index, row in enumerate(rows):
            for column, value in row.items():
                row_indices.append(row_index)
                columns.append(column)
                values.append(value)
        self.a_ub = coo_matrix(
            (values, (row_indices, columns)),
            shape=(len(rows), self.variable_count),
        ).tocsr()
        self.b_ub = np.asarray(rhs)

        equalities: list[np.ndarray] = []
        equality_rhs: list[float] = []
        grand_normalization = np.zeros(self.variable_count)
        grand_normalization[v_index(GRAND)] = 1.0
        equalities.append(grand_normalization)
        equality_rhs.append(1.0)
        x_efficiency = np.zeros(self.variable_count)
        x_efficiency[self.x_start : self.x_start + N] = 1.0
        equalities.append(x_efficiency)
        equality_rhs.append(1.0)
        self.a_eq = np.asarray(equalities)
        self.b_eq = np.asarray(equality_rhs)
        self.bounds = [(None, None)] * GRAND + [(0.0, None)] + [(None, None)] * N

    def objective(
        self, support: tuple[int, ...], weights: tuple[F, ...]
    ) -> np.ndarray:
        objective = np.zeros(self.variable_count)
        for mask, weight_fraction in zip(support, weights, strict=True):
            weight = float(weight_fraction)
            if mask & (mask - 1) == 0 and self.coalition & mask:
                player = mask.bit_length() - 1
                objective[self.x_start + player] -= weight
            else:
                objective[mask - 1] -= weight
                if mask == self.coalition:
                    objective[self.delta_index] -= weight
        # linprog minimizes; subtracting the fixed upper grand worth is done
        # when interpreting the result.  Proper bumps keep that worth at one.
        return objective


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-src",
        type=Path,
        default=Path(
            "/Users/maxf/projects/economics-research/game-theory/"
            "exact-game-monotone-selection/src"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    facets = load_facets(args.project_src)
    vertices = balanced_vertices()
    tested = 0
    global_maximum = -float("inf")
    maximizer: dict[str, object] | None = None
    witness: dict[str, object] | None = None

    # Grand-coalition bumps require a separate varying-efficiency formulation.
    # Every proper coalition is covered here; singleton cases are included.
    for coalition in range(1, GRAND):
        model = PairModel(coalition, facets)
        protected_singletons = {
            1 << player for player in range(N) if coalition >> player & 1
        }
        relevant = [
            vertex
            for vertex in vertices
            if coalition in vertex[0]
            and protected_singletons.intersection(vertex[0])
        ]
        coalition_maximum = -float("inf")
        for support, weights in relevant:
            result = linprog(
                model.objective(support, weights),
                A_ub=model.a_ub,
                b_ub=model.b_ub,
                A_eq=model.a_eq,
                b_eq=model.b_eq,
                bounds=model.bounds,
                method="highs",
            )
            if not result.success:
                raise RuntimeError(result.message)
            violation = -float(result.fun) - 1.0
            tested += 1
            if violation > coalition_maximum:
                coalition_maximum = violation
            if violation > global_maximum:
                global_maximum = violation
                maximizer = {
                    "coalition": coalition,
                    "support": list(support),
                    "weights": [str(weight) for weight in weights],
                    "violation": violation,
                    "delta": float(result.x[model.delta_index]),
                }
            if violation > 1e-8:
                witness = dict(maximizer or {})
                witness["solution"] = [float(value) for value in result.x]
                break
        print(
            json.dumps(
                {
                    "coalition": coalition,
                    "balanced_vertices": len(relevant),
                    "maximum_violation": coalition_maximum,
                }
            ),
            flush=True,
        )
        if witness is not None:
            break

    payload = {
        "status": (
            "nonextendable_exact_pair_found"
            if witness is not None
            else "all_proper_exact_bumps_pointwise_extendable"
        ),
        "facet_count": len(facets),
        "balanced_vertex_count": len(vertices),
        "objectives_tested": tested,
        "maximum_violation": global_maximum,
        "maximizer": maximizer,
        "witness": witness,
        "scope": "five-player monotone normalized exact games; proper bumps",
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if witness is not None else 0


if __name__ == "__main__":
    raise SystemExit(main())
