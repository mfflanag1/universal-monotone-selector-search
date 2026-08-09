#!/usr/bin/env python3
"""Search exact balanced pieces of a mixed terminal objective directly."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog

from n5_mixed_terminal_milp import lower_dual_vertices
from n5_mixed_terminal_search import (
    complete_bipartite_topology,
    complete_bipartite_rectangular_topology,
    GameCone,
    four_cycle_topology,
    lower_expectation_dual,
)
from n5_nonlarge_mixed_cycle import shared_sink_topology


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--weights", type=int, nargs=4, default=(2, 1, 2, 1))
    parser.add_argument(
        "--shape",
        choices=("four-cycle", "shared-sink", "k23", "k33"),
        default="four-cycle",
    )
    parser.add_argument("--k33-masks", type=int, nargs=9)
    parser.add_argument("--k23-masks", type=int, nargs=6)
    parser.add_argument("--samples", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--exhaustive", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    if args.shape == "k33":
        if args.k33_masks is None:
            parser.error("--shape k33 requires --k33-masks")
        topology = complete_bipartite_topology(
            tuple(args.k33_masks), (1,) * 9, 3, 5
        )
    elif args.shape == "k23":
        if args.k23_masks is None:
            parser.error("--shape k23 requires --k23-masks")
        topology = complete_bipartite_rectangular_topology(
            tuple(args.k23_masks), (1,) * 6, 2, 3, 5
        )
    else:
        topology = (
            shared_sink_topology()
            if args.shape == "shared-sink"
            else four_cycle_topology(tuple(args.weights), 5)
        )
    cone = GameCone(5, topology)
    candidate_sets = [
        lower_dual_vertices(divergence, 5)
        for divergence in topology.divergence
    ]
    mixed_nodes = [
        node for node, candidates in enumerate(candidate_sets)
        if len(candidates) > 1
    ]
    if args.exhaustive:
        choices: Any = itertools.product(
            *(range(len(candidate_sets[node])) for node in mixed_nodes)
        )
        sample_count = int(
            np.prod([len(candidate_sets[node]) for node in mixed_nodes])
        )
    else:
        choices = (
            tuple(
                rng.randrange(len(candidate_sets[node]))
                for node in mixed_nodes
            )
            for _ in range(args.samples)
        )
        sample_count = args.samples
    best: dict[str, Any] | None = None
    tested = 0
    for mixed_choice in choices:
        selected = [0] * topology.node_count
        for node, choice in zip(mixed_nodes, mixed_choice, strict=True):
            selected[node] = choice
        objective = np.zeros(cone.variable_count)
        constant = 0.0
        for node, choice in enumerate(selected):
            candidate = candidate_sets[node][choice]
            constant += float(candidate["constant"])
            for coalition, coefficient in candidate["coefficients"].items():
                objective[cone.column(node, coalition)] += float(coefficient)
        result = linprog(
            -objective,
            A_ub=cone.a_ub,
            b_ub=cone.b_ub,
            A_eq=cone.a_eq,
            b_eq=cone.b_eq,
            bounds=[(0.0, 1.0)] * cone.variable_count,
            method="highs",
        )
        tested += 1
        if not result.success:
            continue
        selected_piece_value = constant - float(result.fun)
        if best is not None and selected_piece_value <= float(best["objective"]):
            continue
        games = [
            result.x[node * 32 : (node + 1) * 32]
            for node in range(topology.node_count)
        ]
        exact_values = []
        for game, divergence in zip(
            games, topology.divergence, strict=True
        ):
            dual = lower_expectation_dual(game, divergence, 5)
            if dual is None:
                raise RuntimeError("lower-expectation solve failed")
            exact_values.append(dual[0])
        true_objective = sum(exact_values)
        best = {
            "objective": true_objective,
            "selected_piece_value": selected_piece_value,
            "selected_candidates": selected,
            "node_lower_expectations": exact_values,
            "games": [[float(value) for value in game] for game in games],
        }
        if true_objective > 1e-8:
            break
    payload = {
        "status": (
            "positive_terminal_obstruction_found"
            if best is not None and float(best["objective"]) > 1e-8
            else "no_positive_piece_found"
        ),
        "weights": list(args.weights),
        "shape": args.shape,
        "candidate_counts": [len(candidates) for candidates in candidate_sets],
        "mixed_nodes": mixed_nodes,
        "requested_samples": sample_count,
        "tested_samples": tested,
        "best": best,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
