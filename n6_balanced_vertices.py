#!/usr/bin/env python3
"""Generate all minimal balanced collections on six labelled players."""

from __future__ import annotations

import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np


F = Fraction


def weights_for(support: tuple[int, ...], n: int) -> tuple[F, ...]:
    matrix = np.asarray(
        [[int(mask >> player & 1) for mask in support] for player in range(n)],
        dtype=float,
    )
    if np.linalg.matrix_rank(matrix) != len(support):
        raise ValueError("support is not independent")
    solution, *_ = np.linalg.lstsq(matrix, np.ones(n), rcond=None)
    if np.max(abs(matrix @ solution - 1)) > 1e-8:
        raise ValueError("support is not balanced")
    rational = tuple(F(float(value)).limit_denominator(10000) for value in solution)
    if any(value <= 0 for value in rational):
        raise ValueError("support is not minimally balanced")
    if not all(
        sum(
            rational[index] * int(mask >> player & 1)
            for index, mask in enumerate(support)
        )
        == 1
        for player in range(n)
    ):
        raise ValueError("rational reconstruction failed")
    return rational


def n5_collections() -> dict[tuple[int, ...], tuple[F, ...]]:
    result: dict[tuple[int, ...], tuple[F, ...]] = {}
    incidence = {
        mask: tuple((mask >> player) & 1 for player in range(5))
        for mask in range(1, 32)
    }
    for width in range(1, 6):
        for support in itertools.combinations(range(1, 32), width):
            matrix = np.asarray([incidence[mask] for mask in support], float).T
            if np.linalg.matrix_rank(matrix) != width:
                continue
            solution, *_ = np.linalg.lstsq(matrix, np.ones(5), rcond=None)
            if np.min(solution) <= 1e-10 or np.max(abs(matrix @ solution - 1)) > 1e-9:
                continue
            weights = tuple(F(float(value)).limit_denominator(120) for value in solution)
            if all(
                sum(
                    weights[index] * incidence[mask][player]
                    for index, mask in enumerate(support)
                )
                == 1
                for player in range(5)
            ):
                result[support] = weights
    if len(result) != 1292:
        raise RuntimeError(f"expected 1292 n=5 collections, got {len(result)}")
    return result


def rank(support: tuple[int, ...], n: int = 5) -> int:
    matrix = np.asarray(
        [[(mask >> player) & 1 for mask in support] for player in range(n)],
        dtype=float,
    )
    return int(np.linalg.matrix_rank(matrix))


def generate_n6() -> dict[tuple[int, ...], tuple[F, ...]]:
    old = n5_collections()
    supports: set[tuple[int, ...]] = set()
    new_player = 32

    for support, weights in old.items():
        width = len(support)
        for bits in range(1 << width):
            chosen = [index for index in range(width) if bits >> index & 1]
            subtotal = sum((weights[index] for index in chosen), F(0))
            if subtotal > 1:
                continue
            base_lifted = [
                mask | new_player if bits >> index & 1 else mask
                for index, mask in enumerate(support)
            ]
            lifted = list(base_lifted)
            if subtotal < 1:
                lifted.append(new_player)
            supports.add(tuple(sorted(lifted)))

            remainder = 1 - subtotal
            if remainder <= 0:
                continue
            for delta in range(width):
                if bits >> delta & 1 or weights[delta] <= remainder:
                    continue
                split = list(base_lifted)
                split.append(support[delta] | new_player)
                supports.add(tuple(sorted(split)))

    items = list(old.items())
    for left_index, (left, left_weights) in enumerate(items):
        left_set = set(left)
        for right, right_weights in items[left_index + 1 :]:
            union = tuple(sorted(left_set | set(right)))
            width = len(union)
            if width > 6 or rank(union) != width - 1:
                continue
            left_map = dict(zip(left, left_weights, strict=True))
            right_map = dict(zip(right, right_weights, strict=True))
            mu = tuple(left_map.get(mask, F(0)) for mask in union)
            nu = tuple(right_map.get(mask, F(0)) for mask in union)
            for bits in range(1, (1 << width) - 1):
                mu_sum = sum(
                    (mu[index] for index in range(width) if bits >> index & 1),
                    F(0),
                )
                nu_sum = sum(
                    (nu[index] for index in range(width) if bits >> index & 1),
                    F(0),
                )
                if mu_sum == nu_sum:
                    continue
                parameter = (1 - mu_sum) / (nu_sum - mu_sum)
                if not 0 < parameter < 1:
                    continue
                lifted = tuple(
                    sorted(
                        mask | new_player if bits >> index & 1 else mask
                        for index, mask in enumerate(union)
                    )
                )
                supports.add(lifted)

    result = {}
    invalid = []
    for support in supports:
        try:
            result[support] = weights_for(support, 6)
        except ValueError:
            invalid.append(support)
    if invalid:
        raise RuntimeError(
            f"generated {len(invalid)} invalid supports; first={invalid[0]}"
        )
    if len(result) != 200214:
        raise RuntimeError(
            f"expected 200214 n=6 collections, got {len(result)}"
        )
    return result


def write_json(path: Path) -> None:
    collections = generate_n6()
    payload = [
        [list(support), [str(weight) for weight in weights]]
        for support, weights in sorted(collections.items())
    ]
    path.write_text(json.dumps(payload, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    write_json(args.output)
    print(json.dumps({"status": "complete", "output": str(args.output)}))
