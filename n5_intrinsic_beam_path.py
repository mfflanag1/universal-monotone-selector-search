#!/usr/bin/env python3
"""Find canonical intrinsic exact paths with bounded beam search."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from exactified_obstruction_paths import exact_closure
from family_search import (
    Family,
    box_family_slack_float,
    common_core_gap_float,
    family_slack_float,
    serial_family,
)
from n5_exactified_obstruction_paths import (
    GRAND,
    exact_root,
    null_lift,
)
from n5_facet_search import load_facets


F = Fraction


def beam_path(
    base: tuple[F, ...],
    endpoint: tuple[F, ...],
    differences: tuple[int, ...],
    steps: int,
    cap: F,
    beam_width: int,
    facets: list[tuple[int, ...]],
) -> tuple[list[tuple[int, ...]], dict[tuple[int, ...], F]]:
    origin = (0,) * len(differences)
    destination = (steps,) * len(differences)
    roots = {origin: exact_root(base, facets)}
    if roots[origin] is None or roots[origin] > cap:
        raise RuntimeError("origin exceeds the exact grand cap")
    parents: dict[tuple[int, ...], tuple[int, ...] | None] = {
        origin: None
    }
    frontier = {origin}
    cache = dict(roots)
    for depth in range(1, steps * len(differences) + 1):
        candidates = set()
        predecessors = {}
        for state in frontier:
            for index in range(len(differences)):
                if state[index] == steps:
                    continue
                child = list(state)
                child[index] += 1
                child_state = tuple(child)
                candidates.add(child_state)
                predecessors.setdefault(child_state, state)
        feasible = []
        for state in candidates:
            root = cache.get(state)
            if state not in cache:
                game = list(base)
                for level, coalition in zip(
                    state, differences, strict=True
                ):
                    game[coalition] = (
                        base[coalition]
                        + (endpoint[coalition] - base[coalition])
                        * F(level, steps)
                    )
                root = exact_root(tuple(game), facets)
                cache[state] = root
            if root is None or root > cap:
                continue
            remaining = [steps - level for level in state]
            score = (
                max(remaining),
                sum(value * value for value in remaining),
                root,
                state,
            )
            feasible.append((score, state))
        feasible.sort()
        frontier = {
            state for _, state in feasible[:beam_width]
        }
        for state in frontier:
            roots[state] = cache[state]
            parents[state] = predecessors[state]
        print(
            json.dumps(
                {
                    "depth": depth,
                    "candidates": len(candidates),
                    "feasible": len(feasible),
                    "retained": len(frontier),
                    "destination_retained": destination in frontier,
                }
            ),
            flush=True,
        )
        if destination in frontier:
            break
        if not frontier:
            raise RuntimeError("beam lost every exact path")
    if destination not in parents:
        raise RuntimeError("beam did not reach the endpoint")
    path = []
    state: tuple[int, ...] | None = destination
    while state is not None:
        path.append(state)
        state = parents[state]
    path.reverse()
    return path, roots


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--grand-cap", type=F, default=F(7))
    parser.add_argument("--beam-width", type=int, default=500)
    parser.add_argument("--reuse-first-path", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base4 = tuple(
        map(
            F,
            [0, 0, 0, 0, 0, 3, 3, 1, 0, 3, 3, 1, 3, 1, 5, 6],
        )
    )
    arm_a4 = list(base4)
    arm_a4[7] += 3
    arm_b4 = list(base4)
    arm_b4[14] += 1
    closures = [
        null_lift(exact_closure(game, 4))
        for game in (base4, tuple(arm_a4), tuple(arm_b4))
    ]
    base = closures[0]
    facets = load_facets()[0]
    branch_paths = []
    differences_by_branch = []
    for branch, endpoint in enumerate(closures[1:], start=1):
        differences = tuple(
            coalition
            for coalition in range(1, GRAND)
            if endpoint[coalition] != base[coalition]
        )
        if args.reuse_first_path and branch_paths:
            path = list(branch_paths[0])
            for state in path:
                game = list(base)
                for level, coalition in zip(
                    state, differences, strict=True
                ):
                    game[coalition] = (
                        base[coalition]
                        + (
                            endpoint[coalition]
                            - base[coalition]
                        )
                        * F(level, args.steps)
                    )
                root = exact_root(tuple(game), facets)
                if root is None or root > args.grand_cap:
                    raise RuntimeError(
                        "the first branch path is not symmetric"
                    )
        else:
            path, _ = beam_path(
                base,
                endpoint,
                differences,
                args.steps,
                args.grand_cap,
                args.beam_width,
                facets,
            )
        branch_paths.append(path)
        differences_by_branch.append(differences)
        print(
            json.dumps(
                {
                    "branch": branch,
                    "directions": list(differences),
                    "path_length": len(path),
                }
            ),
            flush=True,
        )

    games = []
    game_index = {}
    edges = []
    edge_set = set()
    depths = []
    for endpoint, differences, path in zip(
        closures[1:],
        differences_by_branch,
        branch_paths,
        strict=True,
    ):
        previous_node = None
        previous_state = None
        for state in path:
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
            exact_game = tuple(game)
            node = game_index.get(exact_game)
            if node is None:
                node = len(games)
                games.append(exact_game)
                depths.append(sum(state))
                game_index[exact_game] = node
            if previous_node is not None and node != previous_node:
                changed = [
                    coalition
                    for coalition, (left, right) in enumerate(
                        zip(
                            games[previous_node],
                            exact_game,
                            strict=True,
                        )
                    )
                    if left != right
                ]
                if len(changed) != 1:
                    raise RuntimeError("beam path step is not one-coordinate")
                edge = (previous_node, node, changed[0])
                if edge not in edge_set:
                    edges.append(edge)
                    edge_set.add(edge)
            previous_node = node
            previous_state = state

    family = Family(
        games,
        edges,
        depths,
        "n5_intrinsic_beam_canonical_paths",
    )
    margin, _ = family_slack_float(games, edges, 5)
    box_margin, _ = box_family_slack_float(games, edges, 5)
    common_gap = common_core_gap_float(games, 5)
    result = {
        "status": "n5_intrinsic_beam_path_complete",
        "n": 5,
        "steps": args.steps,
        "grand_worth": str(args.grand_cap),
        "beam_width": args.beam_width,
        "node_count": len(games),
        "edge_count": len(edges),
        "common_core_budget_gap_float": common_gap,
        "max_min_monotonicity_margin_float": margin,
        "box_max_min_monotonicity_margin_float": box_margin,
        "non_atomic_facet_tax_float": box_margin - margin,
        "branch_paths": [
            [list(state) for state in path] for path in branch_paths
        ],
        "family": serial_family(family, margin, common_gap),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "node_count",
                    "edge_count",
                    "common_core_budget_gap_float",
                    "max_min_monotonicity_margin_float",
                    "box_max_min_monotonicity_margin_float",
                    "non_atomic_facet_tax_float",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
