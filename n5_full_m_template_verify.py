#!/usr/bin/env python3
"""Prove exact feasibility and reachability of the repeated full-m templates."""

from __future__ import annotations

import argparse
import json
from collections import deque
from fractions import Fraction
from pathlib import Path
from typing import Iterable

from n5_facet_search import load_facets
from n5_full_m_recurrence_audit import (
    branch_geometry,
    checkpoint_branch_states,
    load,
    shift,
)


F = Fraction
N = 5
GRAND = 31


def inequality_rows() -> list[tuple[int, ...]]:
    rows = []
    for coalition in range(GRAND + 1):
        for player in range(N):
            if coalition >> player & 1:
                continue
            upper = coalition | (1 << player)
            row = [0] * (GRAND + 1)
            row[upper] = 1
            row[coalition] = -1
            rows.append(tuple(row))
    rows.extend(load_facets()[0])
    return rows


def successors(
    state: tuple[int, ...], allowed: set[tuple[int, ...]]
) -> Iterable[tuple[int, ...]]:
    for coordinate in range(len(state)):
        child = list(state)
        child[coordinate] += 1
        candidate = tuple(child)
        if candidate in allowed:
            yield candidate


def reachable_from(
    seeds: set[tuple[int, ...]], allowed: set[tuple[int, ...]]
) -> set[tuple[int, ...]]:
    seen = set(seeds)
    queue = deque(seeds)
    while queue:
        state = queue.popleft()
        for child in successors(state, allowed):
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return seen


def edge_count(states: set[tuple[int, ...]]) -> int:
    total = 0
    for coordinate in range(10):
        lines: dict[tuple[int, ...], int] = {}
        for state in states:
            line = state[:coordinate] + state[coordinate + 1 :]
            lines[line] = lines.get(line, 0) + 1
        total += sum(length * (length - 1) // 2 for length in lines.values())
    return total


def cross_edge_count(
    left: set[tuple[int, ...]],
    right: set[tuple[int, ...]],
    right_shift: int,
) -> int:
    """Count one-coordinate edges from left to a shifted right template."""
    total = 0
    for coordinate in range(10):
        right_lines: dict[tuple[int, ...], list[int]] = {}
        for state in right:
            shifted_line = tuple(
                value + right_shift
                for index, value in enumerate(state)
                if index != coordinate
            )
            right_lines.setdefault(shifted_line, []).append(
                state[coordinate] + right_shift
            )
        for state in left:
            line = state[:coordinate] + state[coordinate + 1 :]
            total += sum(
                state[coordinate] != right_value
                for right_value in right_lines.get(line, ())
            )
    return total


def parametrically_nonnegative(
    template_name: str,
    slope_m: int,
    slope_shift: int,
    constant: int,
) -> bool:
    if template_name == "P":
        return slope_m >= 0 and 5 * slope_m + constant >= 0
    if template_name == "R":
        return (
            slope_m >= 0
            and slope_m + slope_shift >= 0
            and 6 * slope_m + constant >= 0
        )
    terminal_slope = slope_m + slope_shift
    return terminal_slope >= 0 and 5 * slope_m + constant >= 0


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
    failures = []
    branch_results = []
    total_parametric_checks = 0
    total_common_budget_checks = 0
    origin = (0,) * 10
    common_budget_allocation = (1, 1, 3, 3, 0)
    singleton_witnesses: dict[int, dict[str, object]] = {}
    branch_coordinate_sets_disjoint = not (
        set(differences[0]) & set(differences[1])
    )
    branch_directions_positive = all(
        endpoint[coalition] - base[coalition] > 0
        for endpoint, changed in zip(endpoints, differences, strict=True)
        for coalition in changed
    )
    if not branch_coordinate_sets_disjoint:
        failures.append("the two branches have overlapping proper coordinates")
    if not branch_directions_positive:
        failures.append("a branch direction is not strictly positive")

    for branch in (0, 1):
        states5 = checkpoint_branch_states(
            checkpoints[5][1],
            5,
            branch,
            base,
            endpoints,
            differences,
        )
        states6 = checkpoint_branch_states(
            checkpoints[6][1],
            6,
            branch,
            base,
            endpoints,
            differences,
        )
        early = {state for state in states5 if sum(state) < 24}
        terminal = {state for state in states5 if sum(state) >= 24}
        repeated = {
            state for state in states6 if 24 <= sum(state) < 34
        }
        if any(not 0 <= value <= 5 for state in early for value in state):
            failures.append(f"branch {branch + 1}: P violates level bounds")
        if any(not 0 <= value <= 5 for state in terminal for value in state):
            failures.append(f"branch {branch + 1}: T violates level bounds")
        if any(not 0 <= value <= 6 for state in repeated for value in state):
            failures.append(f"branch {branch + 1}: R violates level bounds")

        endpoint = endpoints[branch]
        changed = differences[branch]
        if any(value.denominator != 1 for value in base):
            raise RuntimeError("base closure is not integral")
        deltas_exact = [
            endpoint[coalition] - base[coalition] for coalition in changed
        ]
        if any(value.denominator != 1 for value in deltas_exact):
            raise RuntimeError("endpoint directions are not integral")
        base_integer = [int(value) for value in base]
        deltas = [int(value) for value in deltas_exact]
        for row_index, row in enumerate(rows):
            slope_m = sum(
                row[coalition] * base_integer[coalition]
                for coalition in range(32)
            )
            direction = [
                row[coalition] * delta
                for coalition, delta in zip(changed, deltas, strict=True)
            ]
            slope_shift = sum(direction)
            grand_offset = 4 * row[GRAND]
            for template_name, states in (
                ("P", early),
                ("R", repeated),
                ("T", terminal),
            ):
                for state in states:
                    constant = grand_offset + sum(
                        coefficient * level
                        for coefficient, level in zip(
                            direction, state, strict=True
                        )
                    )
                    total_parametric_checks += 1
                    valid = parametrically_nonnegative(
                        template_name,
                        slope_m,
                        slope_shift,
                        constant,
                    )
                    if not valid:
                        failures.append(
                            f"branch {branch + 1} {template_name} state "
                            f"{state} fails parametric row {row_index}"
                        )
                        break
                if failures:
                    break
            if failures:
                break

        for coalition in range(1, GRAND):
            budget_value = sum(
                common_budget_allocation[player]
                for player in range(N)
                if coalition >> player & 1
            )
            slope_m = budget_value - base_integer[coalition]
            direction = [
                -delta if changed_coalition == coalition else 0
                for changed_coalition, delta in zip(
                    changed, deltas, strict=True
                )
            ]
            slope_shift = sum(direction)
            for template_name, states in (
                ("P", early),
                ("R", repeated),
                ("T", terminal),
            ):
                for state in states:
                    constant = sum(
                        coefficient * level
                        for coefficient, level in zip(
                            direction, state, strict=True
                        )
                    )
                    total_common_budget_checks += 1
                    if not parametrically_nonnegative(
                        template_name,
                        slope_m,
                        slope_shift,
                        constant,
                    ):
                        failures.append(
                            f"branch {branch + 1} {template_name} state "
                            f"{state} exceeds common budget allocation on "
                            f"coalition {coalition}"
                        )
                        break
                if failures:
                    break
            if failures:
                break

        changed_index = {
            coalition: index for index, coalition in enumerate(changed)
        }
        for player, target in enumerate(common_budget_allocation):
            if player in singleton_witnesses:
                continue
            coalition = 1 << player
            index = changed_index.get(coalition)
            delta = deltas[index] if index is not None else 0
            for template_name, states in (("P", early), ("T", terminal)):
                for state in states:
                    level = state[index] if index is not None else 0
                    if template_name == "P":
                        identity = (
                            target == base_integer[coalition]
                            and delta * level == 0
                        )
                    else:
                        identity = (
                            target
                            == base_integer[coalition] + delta
                            and delta * (level - 5) == 0
                        )
                    if identity:
                        singleton_witnesses[player] = {
                            "branch": branch + 1,
                            "template": template_name,
                            "state": state,
                            "value": target,
                        }
                        break
                if player in singleton_witnesses:
                    break

        p_reachable = reachable_from({origin}, early)
        m5_reachable = reachable_from(p_reachable, early | terminal)
        first_block_reachable = reachable_from(
            p_reachable, early | repeated
        )
        next_block = shift(repeated, 1)
        repeated_reachable = reachable_from(
            repeated, repeated | next_block
        )
        first_terminal = shift(terminal, 1)
        terminal_reachable = reachable_from(
            repeated, repeated | first_terminal
        )
        reachability_checks = {
            "P_from_origin": p_reachable == early,
            "T_from_P_at_m5": terminal <= m5_reachable,
            "first_R_from_P": repeated <= first_block_reachable,
            "next_R_from_previous_R": next_block <= repeated_reachable,
            "shifted_T_from_last_R": first_terminal <= terminal_reachable,
        }
        for name, passed in reachability_checks.items():
            if not passed:
                failures.append(f"branch {branch + 1}: {name} failed")
        family5 = early | terminal
        family6 = early | repeated | shift(terminal, 1)
        family7 = (
            early
            | repeated
            | shift(repeated, 1)
            | shift(terminal, 2)
        )
        edges5 = edge_count(family5)
        edges6 = edge_count(family6)
        edges7 = edge_count(family7)
        edge_increment = edges6 - edges5
        if edges7 - edges6 != edge_increment:
            failures.append(
                f"branch {branch + 1}: repeated edge increment failed"
            )
        if edge_increment != 12766:
            failures.append(
                f"branch {branch + 1}: expected edge increment 12766"
            )
        interaction_offsets = {
            "P_to_R": {
                str(offset): cross_edge_count(early, repeated, offset)
                for offset in range(6)
            },
            "P_to_T": {
                str(offset): cross_edge_count(early, terminal, offset)
                for offset in range(6)
            },
            "R_to_R": {
                str(offset): cross_edge_count(repeated, repeated, offset)
                for offset in range(1, 7)
            },
            "R_to_T": {
                str(offset): cross_edge_count(repeated, terminal, offset)
                for offset in range(1, 7)
            },
        }
        forbidden_nonlocal = {
            "P_to_R": range(1, 6),
            "P_to_T": range(1, 6),
            "R_to_R": range(2, 7),
            "R_to_T": range(2, 7),
        }
        for interaction, offsets in forbidden_nonlocal.items():
            for offset in offsets:
                if interaction_offsets[interaction][str(offset)] != 0:
                    failures.append(
                        f"branch {branch + 1}: nonlocal {interaction} "
                        f"edges occur at offset {offset}"
                    )
        branch_results.append(
            {
                "branch": branch + 1,
                "P_states": len(early),
                "R_states": len(repeated),
                "T_states": len(terminal),
                "edges_at_m5": edges5,
                "edges_at_m6": edges6,
                "edges_at_m7": edges7,
                "edges_per_added_R_block": edge_increment,
                "template_edge_interactions_by_shift": interaction_offsets,
                "nonlocal_edge_interactions_excluded": not any(
                    interaction_offsets[interaction][str(offset)]
                    for interaction, offsets in forbidden_nonlocal.items()
                    for offset in offsets
                ),
                "reachability_checks": reachability_checks,
            }
        )

    if len(singleton_witnesses) != N:
        failures.append("not every singleton maximum has an all-m witness")

    result = {
        "status": "pass" if not failures else "fail",
        "claim": (
            "For every real m>=5, the P, repeated R+k*1, and shifted "
            "T games are monotone exact games; for every integer m>=5 "
            "their unit-coordinate bump graph is reachable from the origin."
        ),
        "inequality_rows_per_game": len(rows),
        "parametric_inequality_checks": total_parametric_checks,
        "parametric_common_budget_checks": total_common_budget_checks,
        "branches": branch_results,
        "all_m_counts": {
            "vertices": "2*(3688 + (m-5)*2604 + 3190) - 1 = 5208m-12285",
            "edges": "2*(32458 + (m-5)*12766) = 25532m-62744",
            "cross_branch_edges": (
                "none beyond the shared origin: the two strictly positive "
                "branch directions have disjoint proper-coalition supports"
            ),
            "branch_coordinate_sets_disjoint": branch_coordinate_sets_disjoint,
            "branch_directions_positive": branch_directions_positive,
        },
        "all_m_global_common_core_gap": {
            "common_budget_threshold": "8",
            "upper_allocation": list(common_budget_allocation),
            "singleton_partition_lower_certificate": {
                str(1 << player): {
                    "weight": "1",
                    "value": value,
                    "witness": singleton_witnesses.get(player),
                }
                for player, value in enumerate(common_budget_allocation)
            },
            "grand_worth": "6+4/m",
            "gap": "2-4/m",
            "proved": len(singleton_witnesses) == N,
        },
        "failures": failures,
        "limitation": (
            "This verifier proves feasibility, reachability, exact affine "
            "counts, and the global gap for the templates. The companion "
            "exclusion verifier proves that no additional state is reachable. "
            "Neither script by itself proves the selector-margin formula."
        ),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
