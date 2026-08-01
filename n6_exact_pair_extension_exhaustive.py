#!/usr/bin/env python3
"""Exhaustively test protected core extension for six-player exact bumps.

All 200,214 minimal balanced collections are reduced by the stabilizer of a
representative bumped coalition of each cardinality.  Each remaining orbit
objective is optimized over the complete monotone six-player exact cone.
"""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from n6_balanced_vertices import weights_for
from n6_exact_pair_extension_sample import GRAND, N, PairModel, load_facets


F = Fraction


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(
        1 << permutation[player]
        for player in range(N)
        if mask >> player & 1
    )


def stabilizer(coalition: int) -> list[tuple[int, ...]]:
    inside = [player for player in range(N) if coalition >> player & 1]
    outside = [player for player in range(N) if not coalition >> player & 1]
    result = []
    for inside_image in itertools.permutations(inside):
        for outside_image in itertools.permutations(outside):
            permutation = list(range(N))
            for source, target in zip(inside, inside_image, strict=True):
                permutation[source] = target
            for source, target in zip(outside, outside_image, strict=True):
                permutation[source] = target
            result.append(tuple(permutation))
    return result


def canonical_support(
    support: tuple[int, ...], permutations: list[tuple[int, ...]]
) -> tuple[int, ...]:
    return min(
        tuple(sorted(permute_mask(mask, permutation) for mask in support))
        for permutation in permutations
    )


def load_collections(path: Path) -> list[tuple[tuple[int, ...], tuple[F, ...]]]:
    raw = json.loads(path.read_text())
    result = [
        (tuple(support), tuple(F(weight) for weight in weights))
        for support, weights in raw
    ]
    if len(result) != 200214:
        raise RuntimeError(f"expected 200214 collections, got {len(result)}")
    return result


def orbit_representatives(
    collections: list[tuple[tuple[int, ...], tuple[F, ...]]], coalition: int
) -> list[tuple[tuple[int, ...], tuple[F, ...]]]:
    protected = {1 << player for player in range(N) if coalition >> player & 1}
    permutations = stabilizer(coalition)
    supports = {
        canonical_support(support, permutations)
        for support, _weights in collections
        if coalition in support and protected.intersection(support)
    }
    return [(support, weights_for(support, N)) for support in sorted(supports)]


def objective(
    model: PairModel,
    support: tuple[int, ...],
    weights: tuple[F, ...],
    singleton_branch: str,
) -> np.ndarray:
    result = np.zeros(model.variable_count)
    for mask, rational_weight in zip(support, weights, strict=True):
        weight = float(rational_weight)
        protected_singleton = mask & (mask - 1) == 0 and mask & model.coalition
        if protected_singleton and singleton_branch == "lower_point":
            result[model.x_start + mask.bit_length() - 1] -= weight
        else:
            result[mask - 1] -= weight
            if mask == model.coalition:
                result[model.delta] -= weight
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collections", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    collections = load_collections(args.collections)
    facets = load_facets(
        Path(
            "/Users/maxf/projects/economics-research/game-theory/"
            "exact-game-monotone-selection/src"
        )
    )
    tested = 0
    labeled_relevant = 0
    maximum = -float("inf")
    maximizer = None
    witness = None
    cardinalities = []

    for size in range(1, N):
        coalition = (1 << size) - 1
        protected = {1 << player for player in range(N) if coalition >> player & 1}
        labeled_count = sum(
            coalition in support and bool(protected.intersection(support))
            for support, _weights in collections
        )
        representatives = orbit_representatives(collections, coalition)
        model = PairModel(coalition, facets)
        branches = ("lower_point", "upper_worth") if size == 1 else ("lower_point",)
        cardinality_maximum = -float("inf")
        for support, weights in representatives:
            for branch in branches:
                result = linprog(
                    objective(model, support, weights, branch),
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
                cardinality_maximum = max(cardinality_maximum, violation)
                if violation > maximum:
                    maximum = violation
                    maximizer = {
                        "coalition": coalition,
                        "coalition_size": size,
                        "branch": branch,
                        "support": list(support),
                        "weights": [str(weight) for weight in weights],
                        "violation": violation,
                        "delta": float(result.x[model.delta]),
                    }
                if violation > 1e-8:
                    witness = dict(maximizer or {})
                    witness["solution"] = [float(value) for value in result.x]
                    break
            if witness is not None:
                break
        cardinalities.append(
            {
                "coalition": coalition,
                "coalition_size": size,
                "labeled_relevant_collections": labeled_count,
                "orbit_representatives": len(representatives),
                "objective_branches": len(branches),
                "maximum_violation": cardinality_maximum,
            }
        )
        labeled_relevant += labeled_count
        print(json.dumps(cardinalities[-1]), flush=True)
        if witness is not None:
            break

    payload = {
        "status": (
            "nonextendable_exact_pair_found"
            if witness is not None
            else "all_proper_exact_bumps_pointwise_extendable"
        ),
        "scope": (
            "six-player monotone normalized exact games; every proper one-coalition "
            "bump; every minimal balanced certificate modulo exact permutation symmetry"
        ),
        "facet_count": len(facets),
        "balanced_collection_count": len(collections),
        "labeled_relevant_collections": labeled_relevant,
        "objectives_tested": tested,
        "maximum_violation": maximum,
        "maximizer": maximizer,
        "witness": witness,
        "cardinalities": cardinalities,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if witness is not None else 0


if __name__ == "__main__":
    raise SystemExit(main())
