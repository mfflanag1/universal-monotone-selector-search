#!/usr/bin/env python3
"""Find one-coordinate bridges between additive-scaled copies of exact games."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from n5_equivariant_boundary_closure import load_quotient
from n5_intrinsic_permutation_union import permute_game


def exact_solution(matrix, target):
    rows = [list(row) + [value] for row, value in zip(matrix, target)]
    pivot_row = 0
    pivots = []
    for column in range(6):
        pivot = next(
            (row for row in range(pivot_row, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            scale = rows[row][column]
            rows[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(
                    rows[row], rows[pivot_row], strict=True
                )
            ]
        pivots.append(column)
        pivot_row += 1
    for row in rows:
        if not any(row[:6]) and row[6]:
            return None
    if len(pivots) != 6:
        return None
    solution = [Fraction(0)] * 6
    for row, column in enumerate(pivots):
        solution[column] = rows[row][6]
    return solution


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--all-permutations", action="store_true")
    parser.add_argument("--nontrivial-only", action="store_true")
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    games, _edges, _stabilizers = load_quotient(args.input)
    permutations = (
        list(itertools.permutations(range(5)))
        if args.all_permutations
        else [tuple(range(5))]
    )
    target_values = np.asarray(
        [[float(value) for value in game] for game in games]
    )
    candidates = []
    best = None
    for right_index, source_right in enumerate(games):
        for permutation in permutations:
            right = permute_game(source_right, permutation)
            right_values = np.asarray([float(value) for value in right])
            for coalition in range(1, 32):
                rows = [mask for mask in range(1, 32) if mask != coalition]
                matrix = np.asarray(
                    [
                        [
                            right_values[mask],
                            *(float(mask >> player & 1) for player in range(5)),
                        ]
                        for mask in rows
                    ]
                )
                solutions, _residuals, rank, _singular = np.linalg.lstsq(
                    matrix, target_values[:, rows].T, rcond=None
                )
                if rank < 6:
                    continue
                residuals = np.max(
                    np.abs(matrix @ solutions - target_values[:, rows].T),
                    axis=0,
                )
                for left_index in np.flatnonzero(residuals < args.tolerance):
                    approx = solutions[:, int(left_index)]
                    delta = (
                        approx[0] * right_values[coalition]
                        + sum(
                            approx[player + 1]
                            for player in range(5)
                            if coalition >> player & 1
                        )
                        - target_values[int(left_index), coalition]
                    )
                    if approx[0] <= args.tolerance or abs(delta) <= args.tolerance:
                        continue
                    exact_matrix = [
                        [
                            right[mask],
                            *(Fraction(mask >> player & 1) for player in range(5)),
                        ]
                        for mask in rows
                    ]
                    solution = exact_solution(
                        exact_matrix,
                        [games[int(left_index)][mask] for mask in rows],
                    )
                    if solution is None or solution[0] <= 0:
                        continue
                    if args.nontrivial_only and solution[0] == 1 and not any(
                        solution[1:]
                    ):
                        continue
                    exact_delta = (
                        solution[0] * right[coalition]
                        + sum(
                            solution[player + 1]
                            for player in range(5)
                            if coalition >> player & 1
                        )
                        - games[int(left_index)][coalition]
                    )
                    if not exact_delta:
                        continue
                    row = {
                        "left": int(left_index),
                        "right": right_index,
                        "permutation": list(permutation),
                        "coalition": coalition,
                        "scale": str(solution[0]),
                        "additive_shift": [str(value) for value in solution[1:]],
                        "delta": str(exact_delta),
                    }
                    candidates.append(row)
                    print(json.dumps(row), flush=True)
                local_best = float(np.min(residuals))
                if best is None or local_best < best[0]:
                    best = (
                        local_best,
                        right_index,
                        permutation,
                        coalition,
                        int(np.argmin(residuals)),
                    )
    output = {
        "status": "affine_bridges_found" if candidates else "no_affine_bridge_found",
        "source": str(args.input),
        "game_count": len(games),
        "permutation_count": len(permutations),
        "candidate_count": len(candidates),
        "candidates": candidates,
        "best_float_residual": None if best is None else best[0],
        "best_float_location": None
        if best is None
        else {
            "right": best[1],
            "permutation": list(best[2]),
            "coalition": best[3],
            "left": best[4],
        },
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: value for key, value in output.items() if key != "candidates"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
