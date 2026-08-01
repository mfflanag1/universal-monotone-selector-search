#!/usr/bin/env python3
"""Prove that the repeated full-m templates contain every reachable state."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from n5_full_m_recurrence_audit import (
    branch_geometry,
    checkpoint_branch_states,
    load,
)
from n5_full_m_template_verify import inequality_rows


F = Fraction
GRAND = 31


def row_coefficients(
    row: tuple[int, ...],
    base: tuple[F, ...],
    changed: list[int],
    endpoint: tuple[F, ...],
) -> tuple[int, int, tuple[int, ...]]:
    base_slope = sum(row[coalition] * base[coalition] for coalition in range(32))
    deltas = tuple(int(endpoint[c] - base[c]) for c in changed)
    direction = tuple(
        row[coalition] * delta
        for coalition, delta in zip(changed, deltas, strict=True)
    )
    if base_slope.denominator != 1:
        raise RuntimeError("nonintegral base row")
    return int(base_slope), 4 * row[GRAND], direction


def numerator(
    coefficients: tuple[int, int, tuple[int, ...]],
    state: tuple[int, ...],
) -> tuple[int, int]:
    slope_m, grand_offset, direction = coefficients
    constant = grand_offset + sum(
        coefficient * level
        for coefficient, level in zip(direction, state, strict=True)
    )
    return slope_m, constant


def negative_on_halfline(slope: int, value_at_start: int) -> bool:
    return slope <= 0 and value_at_start < 0


def witness_fixed(
    state: tuple[int, ...],
    start_m: int,
    coefficients: list[tuple[int, int, tuple[int, ...]]],
) -> int | None:
    for index, coeffs in enumerate(coefficients):
        slope, constant = numerator(coeffs, state)
        if negative_on_halfline(slope, slope * start_m + constant):
            return index
    return None


def witness_shifted(
    state: tuple[int, ...],
    shift_offset: int,
    start_m: int,
    coefficients: list[tuple[int, int, tuple[int, ...]]],
) -> int | None:
    for index, coeffs in enumerate(coefficients):
        slope_m, grand_offset, direction = coeffs
        shift_slope = sum(direction)
        slope = slope_m + shift_slope
        constant = (
            grand_offset
            + sum(
                coefficient * level
                for coefficient, level in zip(direction, state, strict=True)
            )
            + shift_offset * shift_slope
        )
        if negative_on_halfline(slope, slope * start_m + constant):
            return index
    return None


def witness_repeated_interior(
    state: tuple[int, ...],
    coefficients: list[tuple[int, int, tuple[int, ...]]],
) -> int | None:
    """Find a row negative for all m>=7 and 0<=k<=m-7."""
    for index, coeffs in enumerate(coefficients):
        slope_m, grand_offset, direction = coeffs
        slope_k = sum(direction)
        constant = grand_offset + sum(
            coefficient * level
            for coefficient, level in zip(direction, state, strict=True)
        )
        if slope_k >= 0:
            slope = slope_m + slope_k
            value_at_7 = 7 * slope_m + constant
        else:
            slope = slope_m
            value_at_7 = 7 * slope_m + constant
        if negative_on_halfline(slope, value_at_7):
            return index
    return None


def child(state: tuple[int, ...], coordinate: int) -> tuple[int, ...]:
    result = list(state)
    result[coordinate] += 1
    return tuple(result)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m5", type=Path, required=True)
    parser.add_argument("--m6", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    checkpoints = load([args.m5, args.m6])
    if set(checkpoints) != {5, 6}:
        raise ValueError("inputs must be the m=5 and m=6 checkpoints")
    base, endpoints, differences = branch_geometry()
    rows = inequality_rows()
    failures: list[str] = []
    branch_results = []

    for branch in (0, 1):
        states5 = checkpoint_branch_states(
            checkpoints[5][1], 5, branch, base, endpoints, differences
        )
        states6 = checkpoint_branch_states(
            checkpoints[6][1], 6, branch, base, endpoints, differences
        )
        early = {state for state in states5 if sum(state) < 24}
        terminal = {state for state in states5 if sum(state) >= 24}
        repeated = {state for state in states6 if 24 <= sum(state) < 34}
        coefficients = [
            row_coefficients(row, base, differences[branch], endpoints[branch])
            for row in rows
        ]
        counts = {
            "P_internal": 0,
            "P_to_R": 0,
            "P_excluded": 0,
            "R_internal": 0,
            "R_to_next": 0,
            "R_interior_excluded": 0,
            "R_last_internal": 0,
            "R_last_to_T": 0,
            "R_last_excluded": 0,
            "T_internal": 0,
            "T_excluded": 0,
        }
        witnesses: dict[str, dict[str, object]] = {}

        for state in early:
            for coordinate in range(10):
                candidate = child(state, coordinate)
                if candidate in early:
                    counts["P_internal"] += 1
                elif candidate in repeated:
                    counts["P_to_R"] += 1
                else:
                    counts["P_excluded"] += 1
                    row_index = witness_fixed(candidate, 6, coefficients)
                    if row_index is None:
                        failures.append(
                            f"branch {branch + 1}: P child {candidate} has no "
                            "uniform infeasibility witness for m>=6"
                        )
                    else:
                        witnesses.setdefault(
                            "P_excluded",
                            {"state": candidate, "row": row_index},
                        )

        repeated_plus_one = {
            tuple(value + 1 for value in state) for state in repeated
        }
        terminal_plus_one = {
            tuple(value + 1 for value in state) for state in terminal
        }
        for state in repeated:
            for coordinate in range(10):
                candidate = child(state, coordinate)
                if candidate in repeated:
                    counts["R_internal"] += 1
                elif candidate in repeated_plus_one:
                    counts["R_to_next"] += 1
                else:
                    counts["R_interior_excluded"] += 1
                    row_index = witness_repeated_interior(candidate, coefficients)
                    if row_index is None:
                        failures.append(
                            f"branch {branch + 1}: interior R child {candidate} "
                            "has no uniform infeasibility witness"
                        )
                    else:
                        witnesses.setdefault(
                            "R_interior_excluded",
                            {"state": candidate, "row": row_index},
                        )

                if state[coordinate] == 6:
                    continue
                if candidate in repeated:
                    counts["R_last_internal"] += 1
                elif candidate in terminal_plus_one:
                    counts["R_last_to_T"] += 1
                else:
                    counts["R_last_excluded"] += 1
                    row_index = witness_shifted(candidate, -6, 6, coefficients)
                    if row_index is None:
                        failures.append(
                            f"branch {branch + 1}: last R child {candidate} "
                            "has no uniform infeasibility witness"
                        )
                    else:
                        witnesses.setdefault(
                            "R_last_excluded",
                            {"state": candidate, "row": row_index},
                        )

        for state in terminal:
            for coordinate in range(10):
                if state[coordinate] == 5:
                    continue
                candidate = child(state, coordinate)
                if candidate in terminal:
                    counts["T_internal"] += 1
                else:
                    counts["T_excluded"] += 1
                    row_index = witness_shifted(candidate, -5, 5, coefficients)
                    if row_index is None:
                        failures.append(
                            f"branch {branch + 1}: T child {candidate} has no "
                            "uniform infeasibility witness"
                        )
                    else:
                        witnesses.setdefault(
                            "T_excluded",
                            {"state": candidate, "row": row_index},
                        )

        # The two archived base cases ensure equality at m=5 and m=6.
        if len(states5) != len(early | terminal):
            failures.append(f"branch {branch + 1}: m=5 template mismatch")
        expected6 = early | repeated | terminal_plus_one
        if states6 != expected6:
            failures.append(f"branch {branch + 1}: m=6 template mismatch")
        branch_results.append(
            {
                "branch": branch + 1,
                "boundary_counts": counts,
                "sample_infeasibility_witnesses": witnesses,
            }
        )

    result = {
        "status": "pass" if not failures else "fail",
        "claim": (
            "For every integer m>=5, every feasible unit-coordinate child "
            "of a P/R/T template state is again a P/R/T template state. "
            "Together with exact m=5,6 base cases and the reachability "
            "proof, the templates equal the full reachable component."
        ),
        "inequality_rows": len(rows),
        "branches": branch_results,
        "failures": failures,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
