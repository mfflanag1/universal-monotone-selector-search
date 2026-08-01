#!/usr/bin/env python3
"""Sample balanced certificates for protected-extension failure at n=6."""

from __future__ import annotations

import argparse
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


F = Fraction
N = 6
GRAND = 63


def load_facets(project_src: Path):
    sys.path.insert(0, str(project_src))
    from n6_mixed_terminal_search import load_facets as project_load_facets

    return project_load_facets()


def random_balanced_vertex(
    coalition: int, rng: random.Random
) -> tuple[tuple[int, ...], tuple[F, ...]] | None:
    protected = [1 << player for player in range(N) if coalition >> player & 1]
    width = rng.randint(2, N)
    mandatory = {coalition, rng.choice(protected)}
    if len(mandatory) > width:
        return None
    pool = [mask for mask in range(1, GRAND) if mask not in mandatory]
    support = tuple(sorted(mandatory | set(rng.sample(pool, width - len(mandatory)))))
    if len(support) != width:
        return None
    matrix = np.asarray(
        [[(mask >> player) & 1 for mask in support] for player in range(N)],
        dtype=float,
    )
    if np.linalg.matrix_rank(matrix) != width:
        return None
    solution, *_ = np.linalg.lstsq(matrix, np.ones(N), rcond=None)
    if np.min(solution) <= 1e-9 or np.max(np.abs(matrix @ solution - 1.0)) > 1e-8:
        return None
    weights = tuple(F(float(value)).limit_denominator(720) for value in solution)
    if not all(
        sum(
            weights[index] * int(mask >> player & 1)
            for index, mask in enumerate(support)
        )
        == 1
        for player in range(N)
    ):
        return None
    return support, weights


class PairModel:
    def __init__(self, coalition: int, facets) -> None:
        self.coalition = coalition
        self.delta = GRAND
        self.x_start = GRAND + 1
        self.variable_count = GRAND + 1 + N
        rows: list[dict[int, float]] = []
        rhs: list[float] = []

        def add(row: dict[int, float], bound: float = 0.0) -> None:
            rows.append({key: value for key, value in row.items() if value})
            rhs.append(bound)

        for facet in facets:
            lower = {
                mask - 1: -float(coefficient)
                for mask, coefficient in enumerate(facet)
                if mask and coefficient
            }
            add(lower)
            upper = dict(lower)
            if facet[coalition]:
                upper[self.delta] = -float(facet[coalition])
            add(upper)
        for lower_mask in range(GRAND):
            for player in range(N):
                if lower_mask >> player & 1:
                    continue
                upper_mask = lower_mask | (1 << player)
                row = {upper_mask - 1: -1.0}
                if lower_mask:
                    row[lower_mask - 1] = 1.0
                add(row)
                upper_row = dict(row)
                if lower_mask == coalition:
                    upper_row[self.delta] = 1.0
                if upper_mask == coalition:
                    upper_row[self.delta] = -1.0
                add(upper_row)
        for mask in range(1, GRAND):
            row = {mask - 1: 1.0}
            for player in range(N):
                if mask >> player & 1:
                    row[self.x_start + player] = -1.0
            add(row)

        rr, cc, vv = [], [], []
        for row_index, row in enumerate(rows):
            for column, value in row.items():
                rr.append(row_index)
                cc.append(column)
                vv.append(value)
        self.a_ub = coo_matrix(
            (vv, (rr, cc)), shape=(len(rows), self.variable_count)
        ).tocsr()
        self.b_ub = np.asarray(rhs)
        self.a_eq = np.zeros((2, self.variable_count))
        self.a_eq[0, GRAND - 1] = 1.0
        self.a_eq[1, self.x_start : self.x_start + N] = 1.0
        self.b_eq = np.ones(2)
        self.bounds = [(None, None)] * GRAND + [(0.0, None)] + [(None, None)] * N

    def objective(self, support, weights):
        result = np.zeros(self.variable_count)
        for mask, rational_weight in zip(support, weights, strict=True):
            weight = float(rational_weight)
            if mask & (mask - 1) == 0 and mask & self.coalition:
                result[self.x_start + mask.bit_length() - 1] -= weight
            else:
                result[mask - 1] -= weight
                if mask == self.coalition:
                    result[self.delta] -= weight
        return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--project-src",
        type=Path,
        default=Path(
            "/Users/maxf/projects/economics-research/game-theory/"
            "exact-game-monotone-selection/src"
        ),
    )
    args = parser.parse_args()
    facets = load_facets(args.project_src)
    rng = random.Random(args.seed)
    models: dict[int, PairModel] = {}
    seen = set()
    rows = []
    best = None
    attempts = 0
    while len(rows) < args.samples and attempts < args.samples * 1000:
        attempts += 1
        coalition = rng.randrange(1, GRAND)
        candidate = random_balanced_vertex(coalition, rng)
        if candidate is None:
            continue
        support, weights = candidate
        key = (coalition, support, weights)
        if key in seen:
            continue
        seen.add(key)
        model = models.setdefault(coalition, PairModel(coalition, facets))
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
        row = {
            "coalition": coalition,
            "support": list(support),
            "weights": [str(weight) for weight in weights],
            "violation": violation,
            "delta": float(result.x[model.delta]),
        }
        rows.append(row)
        if best is None or violation > best["violation"]:
            best = row | {"solution": [float(value) for value in result.x]}
            print(json.dumps({"tested": len(rows), "best": row}), flush=True)
        if violation > 1e-8:
            break
    payload = {
        "status": (
            "nonextendable_exact_pair_found"
            if best is not None and best["violation"] > 1e-8
            else "none_found"
        ),
        "facet_count": len(facets),
        "samples": len(rows),
        "attempts": attempts,
        "best": best,
        "rows": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "samples", "best")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
