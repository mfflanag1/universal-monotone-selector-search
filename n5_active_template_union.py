#!/usr/bin/env python3
"""Lift the union of archived active P/R/T states to a target resolution."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from family_search import Family, common_core_gap_float, serial_family
from n5_active_component_recurrence_audit import infer_m, state_for_game
from n5_full_m_recurrence_audit import (
    branch_geometry,
    checkpoint_branch_states,
    load,
)
from n5_m5_bridge_closure import induced_edges
from n5_sparse_family_lp import sparse_family_slack_float


F = Fraction
GRAND = 31


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m5", type=Path, required=True)
    parser.add_argument("--m6", type=Path, required=True)
    parser.add_argument("--archive", type=Path, action="append", required=True)
    parser.add_argument("--target-m", type=int, required=True)
    parser.add_argument(
        "--method",
        choices=("highs", "highs-ds", "highs-ipm"),
        default="highs-ipm",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.target_m < 5:
        raise ValueError("target m must be at least 5")

    checkpoints = load([args.m5, args.m6])
    base, endpoints, differences = branch_geometry()
    templates = []
    selected = []
    for branch in (0, 1):
        states5 = checkpoint_branch_states(
            checkpoints[5][1], 5, branch, base, endpoints, differences
        )
        states6 = checkpoint_branch_states(
            checkpoints[6][1], 6, branch, base, endpoints, differences
        )
        templates.append(
            {
                "P": {state for state in states5 if sum(state) < 24},
                "R": {state for state in states6 if 24 <= sum(state) < 34},
                "T": {state for state in states5 if sum(state) >= 24},
            }
        )
        selected.append({"P": set(), "R": set(), "T": set()})

    failures = []
    for path in args.archive:
        source_m = infer_m(path)
        payload = json.loads(path.read_text())
        for raw_game in payload["family"]["games"]:
            game = tuple(F(value) for value in raw_game)
            branch, state = state_for_game(
                game, source_m, base, endpoints, differences
            )
            if branch == 0:
                continue
            template = templates[branch - 1]
            if state in template["P"]:
                selected[branch - 1]["P"].add(state)
                continue
            terminal = tuple(value - (source_m - 5) for value in state)
            if terminal in template["T"]:
                selected[branch - 1]["T"].add(terminal)
                continue
            for block in range(source_m - 5):
                repeated = tuple(value - block for value in state)
                if repeated in template["R"]:
                    selected[branch - 1]["R"].add(repeated)
                    break
            else:
                failures.append(
                    f"unclassified source state m={source_m} branch={branch} {state}"
                )
    if failures:
        raise RuntimeError(failures[0])

    target_m = args.target_m
    grand_worth = F(6) + F(4, target_m)
    games = []
    category_counts = []
    for branch in (0, 1):
        states = set(selected[branch]["P"])
        for block in range(target_m - 5):
            states.update(
                tuple(value + block for value in state)
                for state in selected[branch]["R"]
            )
        states.update(
            tuple(value + target_m - 5 for value in state)
            for state in selected[branch]["T"]
        )
        endpoint = endpoints[branch]
        for state in states:
            game = list(base)
            for level, coalition in zip(
                state, differences[branch], strict=True
            ):
                game[coalition] += (
                    endpoint[coalition] - base[coalition]
                ) * F(level, target_m)
            game[GRAND] = grand_worth
            games.append(tuple(game))
        category_counts.append(
            {
                "branch": branch + 1,
                "P_base_states": len(selected[branch]["P"]),
                "R_base_states": len(selected[branch]["R"]),
                "T_base_states": len(selected[branch]["T"]),
                "instantiated_states": len(states),
            }
        )
    games = list(dict.fromkeys(games))
    edges = induced_edges(games, GRAND)
    margin = sparse_family_slack_float(games, edges, 5, args.method)
    if margin is None:
        raise RuntimeError("target union LP failed")
    gap = common_core_gap_float(games, 5)
    family = Family(
        games,
        edges,
        [0] * len(games),
        "n5_active_template_union",
    )
    result = {
        "status": "n5_active_template_union_complete",
        "n": 5,
        "target_m": target_m,
        "grand_worth": str(grand_worth),
        "source_archives": [str(path) for path in args.archive],
        "category_counts": category_counts,
        "node_count": len(games),
        "edge_count": len(edges),
        "common_core_budget_gap_float": gap,
        "max_min_monotonicity_margin_float": margin,
        "predicted_full_margin": str(
            F(9 * target_m - 28, 4 * target_m * (16 * target_m - 49))
        ),
        "family": serial_family(family, margin, gap),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {key: result[key] for key in result if key != "family"},
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
