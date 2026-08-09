#!/usr/bin/env python3
"""Exhaust the symmetry copies of sharp local-premium branches on a topology."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

from n5_mixed_terminal_search import GameCone, Topology, lower_expectation_dual


F = Fraction
N = 5


def branch_copies(
    divergence: tuple[int, ...], archives: list[dict[str, Any]]
) -> list[tuple[F, tuple[tuple[int, F], ...]]]:
    actual_levels = sorted(set(divergence))
    actual_counts = tuple(divergence.count(level) for level in actual_levels)
    match = next(
        result
        for archive in archives
        for result in archive["results"]
        if tuple(result["level_multiplicities"]) == actual_counts
    )
    canonical = tuple(F(value) for value in match["divergence"])
    canonical_levels = sorted(set(canonical))
    if len(canonical_levels) != len(actual_levels):
        raise RuntimeError("premium archive has the wrong level count")
    scale = F(actual_levels[-1] - actual_levels[0]) / (
        canonical_levels[-1] - canonical_levels[0]
    )
    translation = F(actual_levels[0]) - scale * canonical_levels[0]
    entry = match["certificates"][int(match["witness_branch"])]
    canonical_beta = F(entry["beta"])
    canonical_terms = tuple(
        (int(coalition), F(weight))
        for coalition, weight in entry["terms"]
    )
    canonical_groups = [
        [player for player, value in enumerate(canonical) if value == level]
        for level in canonical_levels
    ]
    actual_groups = [
        [player for player, value in enumerate(divergence) if value == level]
        for level in actual_levels
    ]
    copies = set()
    for targets in itertools.product(
        *(itertools.permutations(group) for group in actual_groups)
    ):
        permutation = {}
        for canonical_group, target_group in zip(
            canonical_groups, targets, strict=True
        ):
            permutation.update(zip(canonical_group, target_group, strict=True))

        def permute_mask(mask: int) -> int:
            return sum(
                1 << permutation[player]
                for player in range(N)
                if mask >> player & 1
            )

        copies.add(
            (
                scale * canonical_beta + translation,
                tuple(
                    sorted(
                        (permute_mask(coalition), scale * weight)
                        for coalition, weight in canonical_terms
                    )
                ),
            )
        )
    return sorted(copies, key=str)


def topology_from_report(report: dict[str, Any]) -> Topology:
    terminals = report["terminals"]
    node_count = max(int(row["node"]) for row in terminals) + 1
    divergences = [[0] * N for _ in range(node_count)]
    for row in terminals:
        divergences[int(row["node"])] = [
            int(value) for value in row["scaled_divergence"]
        ]
    masks: dict[tuple[int, int], int] = {}
    for path in report["terminal_path_decomposition"]:
        key = int(path["source"]), int(path["sink"])
        masks[key] = masks.get(key, 0) | 1 << int(path["player"])
    return Topology(
        node_count,
        tuple(
            (source, sink, mask)
            for (source, sink), mask in sorted(masks.items())
        ),
        tuple(tuple(vector) for vector in divergences),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("--three-level", type=Path, required=True)
    parser.add_argument("--four-level", type=Path, required=True)
    parser.add_argument("--variable-grand", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(args.terminal_report.read_text())
    archives = [
        json.loads(args.three_level.read_text()),
        json.loads(args.four_level.read_text()),
    ]
    topology = topology_from_report(report)
    copies = [branch_copies(divergence, archives) for divergence in topology.divergence]
    combination_count = int(np.prod([len(rows) for rows in copies]))
    cone = GameCone(N, topology, variable_grand=args.variable_grand)
    best = None
    tested = 0
    for selected in itertools.product(*copies):
        objective = np.zeros(cone.variable_count)
        for node, (beta, terms) in enumerate(selected):
            objective[cone.column(node, (1 << N) - 1)] += float(beta)
            for coalition, weight in terms:
                objective[cone.column(node, coalition)] += float(weight)
        point = cone.maximize(objective)
        tested += 1
        if point is None:
            continue
        piece_value = float(objective @ point)
        true_values = []
        for node, divergence in enumerate(topology.divergence):
            game = point[node * 32 : (node + 1) * 32]
            solved = lower_expectation_dual(game, divergence, N)
            if solved is None:
                raise RuntimeError("lower-expectation solve failed")
            true_values.append(solved[0])
        true_value = sum(true_values)
        if best is None or true_value > float(best["true_objective"]):
            best = {
                "true_objective": true_value,
                "selected_piece_objective": piece_value,
                "node_lower_expectations": true_values,
                "selected_representations": [
                    {
                        "beta": str(beta),
                        "terms": [[coalition, str(weight)] for coalition, weight in terms],
                    }
                    for beta, terms in selected
                ],
                "games": [
                    [float(value) for value in point[node * 32 : (node + 1) * 32]]
                    for node in range(topology.node_count)
                ],
            }
        if true_value > 1e-8:
            break
    payload = {
        "status": (
            "positive_terminal_obstruction_found"
            if best is not None and float(best["true_objective"]) > 1e-8
            else "no_positive_aligned_premium_piece_found"
        ),
        "source": str(args.terminal_report),
        "variable_grand": args.variable_grand,
        "copy_counts": [len(rows) for rows in copies],
        "combination_count": combination_count,
        "tested": tested,
        "best": best,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: value
                for key, value in payload.items()
                if key != "best"
            }
            | {
                "best_true_objective": None
                if best is None
                else best["true_objective"],
                "best_piece_objective": None
                if best is None
                else best["selected_piece_objective"],
            },
            indent=2,
        )
    )
    return 1 if payload["status"] == "positive_terminal_obstruction_found" else 0


if __name__ == "__main__":
    raise SystemExit(main())
