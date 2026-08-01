#!/usr/bin/env python3
"""Audit the observed full-grid recurrence and its exact margin laws."""

from __future__ import annotations

import argparse
import json
import re
from fractions import Fraction
from pathlib import Path
from typing import Any

from exactified_obstruction_paths import exact_closure
from n5_exactified_obstruction_paths import null_lift


F = Fraction
REPEATED_LAYER_BLOCK = [305, 304, 299, 265, 247, 226, 224, 220, 241, 273]


def resolution(path: Path, payload: dict[str, Any]) -> int:
    if "steps" in payload:
        return int(payload["steps"])
    match = re.search(r"grid(\d+)", path.name)
    if match is None:
        raise ValueError(f"cannot infer resolution from {path}")
    return int(match.group(1))


def load(paths: list[Path]) -> dict[int, tuple[Path, dict[str, Any]]]:
    result = {}
    for path in paths:
        payload = json.loads(path.read_text())
        result[resolution(path, payload)] = (path, payload)
    return result


def branch_geometry() -> tuple[
    tuple[F, ...], list[tuple[F, ...]], list[list[int]]
]:
    base4 = tuple(
        map(
            F,
            [0, 0, 0, 0, 0, 3, 3, 1, 0, 3, 3, 1, 3, 1, 5, 6],
        )
    )
    endpoints4 = [list(base4), list(base4)]
    endpoints4[0][7] += 3
    endpoints4[1][14] += 1
    closures = [
        null_lift(exact_closure(game, 4))
        for game in (
            base4,
            tuple(endpoints4[0]),
            tuple(endpoints4[1]),
        )
    ]
    differences = [
        [
            coalition
            for coalition in range(1, 31)
            if closures[branch][coalition] != closures[0][coalition]
        ]
        for branch in (1, 2)
    ]
    return closures[0], closures[1:], differences


def checkpoint_branch_states(
    payload: dict[str, Any],
    m: int,
    branch: int,
    base: tuple[F, ...],
    endpoints: list[tuple[F, ...]],
    differences: list[list[int]],
) -> set[tuple[int, ...]]:
    result = set()
    own = differences[branch]
    other = differences[1 - branch]
    endpoint = endpoints[branch]
    for raw_game in payload["family"]["games"]:
        game = tuple(F(value) for value in raw_game)
        state = tuple(
            int(
                (game[coalition] - base[coalition])
                * m
                / (endpoint[coalition] - base[coalition])
            )
            for coalition in own
        )
        if not all(
            base[coalition]
            + (endpoint[coalition] - base[coalition]) * F(level, m)
            == game[coalition]
            for level, coalition in zip(state, own, strict=True)
        ):
            continue
        if any(game[coalition] != base[coalition] for coalition in other):
            continue
        result.add(state)
    return result


def shift(
    states: set[tuple[int, ...]], amount: int
) -> set[tuple[int, ...]]:
    return {
        tuple(level + amount for level in state)
        for state in states
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, action="append", default=[])
    parser.add_argument("--full-float", type=Path, action="append", default=[])
    parser.add_argument("--active-exact", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    checkpoints = load(args.checkpoint)
    full_float = load(args.full_float)
    active_exact = load(args.active_exact)
    failures: list[str] = []
    checkpoint_rows = []
    state_recurrence_rows = []
    prefix = suffix = None
    if checkpoints:
        base_m = min(checkpoints)
        if base_m != 5:
            failures.append("layer recurrence audit requires the m=5 checkpoint")
        else:
            base_layers = checkpoints[5][1]["layer_counts"]
            prefix = base_layers[:24]
            suffix = base_layers[24:]
            if len(suffix) != 27:
                failures.append("m=5 layer suffix does not have length 27")
        for m, (path, payload) in sorted(checkpoints.items()):
            expected_layers = (
                prefix
                + REPEATED_LAYER_BLOCK * (m - 5)
                + suffix
                if prefix is not None and suffix is not None and m >= 5
                else None
            )
            layers_match = payload["layer_counts"] == expected_layers
            reachable_expected = 2604 * m - 6142
            nodes_expected = 5208 * m - 12285
            edges_expected = 25532 * m - 62744
            row = {
                "m": m,
                "source": str(path),
                "layers_match": layers_match,
                "reachable": int(payload["reachable_states_per_branch"]),
                "reachable_expected": reachable_expected,
                "nodes": int(payload["node_count"]),
                "nodes_expected": nodes_expected,
                "edges": int(payload["edge_count"]),
                "edges_expected": edges_expected,
            }
            checkpoint_rows.append(row)
            if not layers_match:
                failures.append(f"m={m} layer histogram breaks the recurrence")
            for key in ("reachable", "nodes", "edges"):
                if row[key] != row[f"{key}_expected"]:
                    failures.append(f"m={m} {key} breaks the affine law")

        if 5 in checkpoints and 6 in checkpoints:
            base, endpoints, differences = branch_geometry()
            states5 = [
                checkpoint_branch_states(
                    checkpoints[5][1],
                    5,
                    branch,
                    base,
                    endpoints,
                    differences,
                )
                for branch in (0, 1)
            ]
            states6 = [
                checkpoint_branch_states(
                    checkpoints[6][1],
                    6,
                    branch,
                    base,
                    endpoints,
                    differences,
                )
                for branch in (0, 1)
            ]
            templates = []
            for branch in (0, 1):
                early = {state for state in states5[branch] if sum(state) < 24}
                terminal = {
                    state for state in states5[branch] if sum(state) >= 24
                }
                repeated = {
                    state
                    for state in states6[branch]
                    if 24 <= sum(state) < 34
                }
                templates.append((early, repeated, terminal))
            for m, (path, payload) in sorted(checkpoints.items()):
                branch_rows = []
                for branch in (0, 1):
                    observed = checkpoint_branch_states(
                        payload,
                        m,
                        branch,
                        base,
                        endpoints,
                        differences,
                    )
                    early, repeated, terminal = templates[branch]
                    expected = set(early) | shift(terminal, m - 5)
                    for copy in range(m - 5):
                        expected |= shift(repeated, copy)
                    matches = observed == expected
                    branch_rows.append(
                        {
                            "branch": branch + 1,
                            "observed_states": len(observed),
                            "expected_states": len(expected),
                            "matches": matches,
                        }
                    )
                    if not matches:
                        failures.append(
                            f"m={m} branch {branch + 1} breaks the "
                            "P + translated-R + translated-T state recurrence"
                        )
                state_recurrence_rows.append(
                    {
                        "m": m,
                        "source": str(path),
                        "branches": branch_rows,
                    }
                )
        elif checkpoints:
            failures.append("state recurrence audit requires m=5 and m=6")

    full_rows = []
    for m, (path, payload) in sorted(full_float.items()):
        observed = F(payload["common_core_budget_gap_float"]).limit_denominator(
            1_000_000
        )
        expected = F(2) - F(4, m)
        full_rows.append(
            {
                "m": m,
                "source": str(path),
                "global_gap": str(observed),
                "global_gap_expected": str(expected),
                "matches": observed == expected,
            }
        )
        if observed != expected:
            failures.append(f"m={m} global gap breaks 2-4/m")

    active_rows = []
    for m, (path, payload) in sorted(active_exact.items()):
        margin = F(payload["max_min_monotonicity_margin_exact"])
        margin_expected = F(9 * m - 28, 4 * m * (16 * m - 49))
        gap = F(payload["common_core_budget_gap_exact"])
        gap_expected = F(2) - F(12, m)
        active_rows.append(
            {
                "m": m,
                "source": str(path),
                "component_gap": str(gap),
                "component_gap_expected": str(gap_expected),
                "margin": str(margin),
                "margin_expected": str(margin_expected),
                "gap_matches": gap == gap_expected,
                "margin_matches": margin == margin_expected,
            }
        )
        if gap != gap_expected:
            failures.append(f"m={m} component gap breaks 2-12/m")
        if margin != margin_expected:
            failures.append(f"m={m} margin breaks the rational law")

    tested = sorted(set(checkpoints) | set(full_float) | set(active_exact))
    next_m = max(tested) + 1 if tested else 10
    result = {
        "status": "pass" if not failures else "fail",
        "tested_resolutions": tested,
        "repeated_layer_block": REPEATED_LAYER_BLOCK,
        "repeated_layer_block_sum": sum(REPEATED_LAYER_BLOCK),
        "checkpoint_rows": checkpoint_rows,
        "state_recurrence": {
            "early_template_states_per_branch": 3688,
            "repeated_template_states_per_branch": 2604,
            "terminal_template_states_per_branch": 3190,
            "rows": state_recurrence_rows,
        },
        "full_family_rows": full_rows,
        "active_component_rows": active_rows,
        "next_prediction": {
            "m": next_m,
            "reachable_states_per_branch": 2604 * next_m - 6142,
            "full_nodes": 5208 * next_m - 12285,
            "full_edges": 25532 * next_m - 62744,
            "global_gap": str(F(2) - F(4, next_m)),
            "component_gap": str(F(2) - F(12, next_m)),
            "margin": str(
                F(
                    9 * next_m - 28,
                    4 * next_m * (16 * next_m - 49),
                )
            ),
        },
        "failures": failures,
        "interpretation": (
            "This is the exact finite archive audit. The companion template "
            "and exclusion verifiers supply the induction proof for all "
            "integer m>=5."
        ),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
