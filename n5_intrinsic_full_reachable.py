#!/usr/bin/env python3
"""Solve the full intrinsic exact component at a fixed grand cap."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from exactified_obstruction_paths import exact_closure
from family_search import (
    Family,
    common_core_gap_float,
    serial_family,
)
from n5_exactified_obstruction_paths import (
    GRAND,
    exact_root,
    null_lift,
)
from n5_facet_search import load_facets
from n5_m5_bridge_closure import induced_edges
from n5_sparse_dual_support import dual_envelope_support
from n5_sparse_family_lp import sparse_family_result


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--grand-cap", type=F, required=True)
    parser.add_argument(
        "--method",
        choices=("highs", "highs-ds", "highs-ipm"),
        default="highs",
    )
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--support-output", type=Path)
    parser.add_argument("--support-threshold", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

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
    base = closures[0]
    facets = load_facets()[0]
    differences_by_branch = [
        tuple(
            coalition
            for coalition in range(1, GRAND)
            if endpoint[coalition] != base[coalition]
        )
        for endpoint in closures[1:]
    ]

    def game_for_state(
        endpoint: tuple[F, ...],
        differences: tuple[int, ...],
        state: tuple[int, ...],
    ) -> tuple[F, ...]:
        game = list(base)
        for level, coalition in zip(
            state, differences, strict=True
        ):
            game[coalition] = (
                base[coalition]
                + (endpoint[coalition] - base[coalition])
                * F(level, args.steps)
            )
        game[GRAND] = args.grand_cap
        return tuple(game)

    differences = differences_by_branch[0]
    origin = (0,) * len(differences)
    destination = (args.steps,) * len(differences)
    frontier = {origin}
    reachable = {origin}
    layer_counts = [1]
    for depth in range(1, args.steps * len(differences) + 1):
        candidates = {
            tuple(
                level + (index == direction)
                for index, level in enumerate(state)
            )
            for state in frontier
            for direction in range(len(differences))
            if state[direction] < args.steps
        }
        feasible = set()
        for state in candidates:
            root = exact_root(
                game_for_state(
                    closures[1], differences, state
                ),
                facets,
            )
            if root is not None and root <= args.grand_cap:
                feasible.add(state)
        frontier = feasible
        reachable.update(frontier)
        layer_counts.append(len(frontier))
        print(
            json.dumps(
                {
                    "depth": depth,
                    "candidates": len(candidates),
                    "feasible": len(frontier),
                    "reachable": len(reachable),
                    "destination_reached": destination in frontier,
                }
            ),
            flush=True,
        )
        if not frontier:
            raise RuntimeError("reachable component died")
    if destination not in reachable:
        raise RuntimeError("endpoint is not reachable")

    games = []
    for endpoint, branch_differences in zip(
        closures[1:], differences_by_branch, strict=True
    ):
        games.extend(
            game_for_state(endpoint, branch_differences, state)
            for state in sorted(reachable, key=lambda value: (sum(value), value))
        )
    games = list(dict.fromkeys(games))
    edges = induced_edges(games, GRAND)
    family = Family(
        games,
        edges,
        [0] * len(games),
        "n5_intrinsic_full_reachable_component",
    )
    print(
        json.dumps(
            {
                "stage": "induced_family_complete",
                "node_count": len(games),
                "edge_count": len(edges),
            }
        ),
        flush=True,
    )
    if args.checkpoint is not None:
        checkpoint = {
            "status": "n5_intrinsic_full_reachable_checkpoint",
            "n": 5,
            "steps": args.steps,
            "grand_worth": str(args.grand_cap),
            "reachable_states_per_branch": len(reachable),
            "layer_counts": layer_counts,
            "node_count": len(games),
            "edge_count": len(edges),
            "family": serial_family(family, None, None),
        }
        args.checkpoint.write_text(
            json.dumps(checkpoint, indent=2) + "\n"
        )
        print(
            json.dumps(
                {
                    "stage": "checkpoint_written",
                    "path": str(args.checkpoint),
                }
            ),
            flush=True,
        )
    solve_result = sparse_family_result(
        games, edges, 5, method=args.method
    )
    if not solve_result.success:
        raise RuntimeError(solve_result.message)
    margin = float(solve_result.x[len(games) * 5])
    common_gap = common_core_gap_float(games, 5)
    if args.support_output is not None:
        support = dual_envelope_support(
            games,
            edges,
            5,
            solve_result,
            str(args.output),
            args.support_threshold,
            args.method,
        )
        args.support_output.write_text(
            json.dumps(support, indent=2) + "\n"
        )
        print(
            json.dumps(
                {
                    "stage": "support_written",
                    "path": str(args.support_output),
                    "node_count": support["node_count"],
                    "edge_count": support["edge_count"],
                    "max_min_monotonicity_margin_float": support[
                        "max_min_monotonicity_margin_float"
                    ],
                }
            ),
            flush=True,
        )
    result = {
        "status": "n5_intrinsic_full_reachable_complete",
        "n": 5,
        "steps": args.steps,
        "grand_worth": str(args.grand_cap),
        "method": args.method,
        "reachable_states_per_branch": len(reachable),
        "layer_counts": layer_counts,
        "node_count": len(games),
        "edge_count": len(edges),
        "common_core_budget_gap_float": common_gap,
        "max_min_monotonicity_margin_float": margin,
        "family": serial_family(family, margin, common_gap),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "reachable_states_per_branch",
                    "node_count",
                    "edge_count",
                    "common_core_budget_gap_float",
                    "max_min_monotonicity_margin_float",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
