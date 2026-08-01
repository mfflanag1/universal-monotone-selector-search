#!/usr/bin/env python3
"""Search the sharpest cycle-rank-two K2,3 protected block flows."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from pathlib import Path
from typing import Any

from n5_mixed_terminal_search import (
    complete_bipartite_rectangular_topology,
    optimize_topology,
)


N = 5
GRAND = (1 << N) - 1
MAXIMUM_PREMIUM_TYPES = {(1, 2, 2), (1, 3, 1), (2, 2, 1)}


def level_multiplicities(
    left: int, right: int, *, negative: bool = False
) -> tuple[int, int, int]:
    intersection = (left & right).bit_count()
    symmetric_difference = (left ^ right).bit_count()
    neither = N - (left | right).bit_count()
    return (
        (intersection, symmetric_difference, neither)
        if negative
        else (neither, symmetric_difference, intersection)
    )


def is_sharp_source_row(row: tuple[int, ...]) -> bool:
    incidence = [
        sum(int(mask >> player & 1) for mask in row)
        for player in range(N)
    ]
    return len(set(incidence)) == 4 and all(
        level_multiplicities(row[left], row[right])
        in MAXIMUM_PREMIUM_TYPES
        for left, right in itertools.combinations(range(3), 2)
    )


def is_all_pairs_sharp(masks: tuple[int, ...]) -> bool:
    if not all(is_sharp_source_row(row) for row in (masks[:3], masks[3:])):
        return False
    return all(
        level_multiplicities(masks[sink], masks[3 + sink], negative=True)
        in MAXIMUM_PREMIUM_TYPES
        for sink in range(3)
    )


def draw_candidate(rng: random.Random) -> tuple[int, ...]:
    while True:
        masks = tuple(rng.sample(range(1, GRAND), 6))
        if not is_all_pairs_sharp(masks):
            continue
        if masks[0] & masks[1] & masks[2] & masks[3] & masks[4] & masks[5]:
            continue
        if masks[0] | masks[1] | masks[2] | masks[3] | masks[4] | masks[5] != GRAND:
            continue
        return masks


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(
        1 << permutation[player]
        for player in range(N)
        if mask >> player & 1
    )


def canonical_topology(masks: tuple[int, ...]) -> tuple[int, ...]:
    rows = masks[:3], masks[3:]
    best: tuple[int, ...] | None = None
    for player_permutation in itertools.permutations(range(N)):
        permuted = tuple(
            tuple(permute_mask(mask, player_permutation) for mask in row)
            for row in rows
        )
        for swap_sources in (False, True):
            ordered_rows = (
                (permuted[1], permuted[0])
                if swap_sources
                else permuted
            )
            for sink_permutation in itertools.permutations(range(3)):
                candidate = tuple(
                    ordered_rows[0][sink] for sink in sink_permutation
                ) + tuple(
                    ordered_rows[1][sink] for sink in sink_permutation
                )
                if best is None or candidate < best:
                    best = candidate
    assert best is not None
    return best


def source_disagreement_coordinates(
    masks: tuple[int, ...]
) -> list[dict[str, Any]]:
    first, second = masks[:3], masks[3:]
    coordinates: list[dict[str, Any]] = []
    for coalition in range(1, GRAND):
        directions: set[int] = set()
        possible = True
        for left, right in zip(first, second, strict=True):
            contains_left = coalition & left == left
            contains_right = coalition & right == right
            if not contains_left and not contains_right:
                possible = False
                break
            if contains_left and not contains_right:
                directions.add(1)
            elif contains_right and not contains_left:
                directions.add(-1)
        if possible and len(directions) <= 1:
            coordinates.append(
                {
                    "coalition": coalition,
                    "relation": (
                        "source_0_le_source_1"
                        if directions == {1}
                        else "source_1_le_source_0"
                        if directions == {-1}
                        else "unconstrained"
                    ),
                }
            )
    return coordinates


def exhaustive_orbits() -> tuple[list[tuple[int, ...]], int]:
    rows = [
        row
        for row in itertools.combinations(range(1, GRAND), 3)
        if is_sharp_source_row(row)
    ]
    representatives: set[tuple[int, ...]] = set()
    sink_normalized_count = 0
    for first in rows:
        for second_set in rows:
            for second in itertools.permutations(second_set):
                masks = first + second
                if not all(
                    level_multiplicities(
                        masks[sink], masks[3 + sink], negative=True
                    )
                    in MAXIMUM_PREMIUM_TYPES
                    for sink in range(3)
                ):
                    continue
                if not is_all_pairs_sharp(masks):
                    continue
                if masks[0] & masks[1] & masks[2] & masks[3] & masks[4] & masks[5]:
                    continue
                if masks[0] | masks[1] | masks[2] | masks[3] | masks[4] | masks[5] != GRAND:
                    continue
                sink_normalized_count += 1
                representatives.add(canonical_topology(masks))
    ordered = sorted(representatives)
    for masks in ordered:
        disagreement = source_disagreement_coordinates(masks)
        assert len(disagreement) == 2
        assert all(row["coalition"].bit_count() == N - 1 for row in disagreement)
        assert {row["relation"] for row in disagreement} == {
            "source_0_le_source_1",
            "source_1_le_source_0",
        }
    return ordered, sink_normalized_count


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--starts", type=int, default=4)
    parser.add_argument("--iterations", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260821)
    parser.add_argument("--variable-grand", action="store_true")
    parser.add_argument("--exhaustive-orbits", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    if args.exhaustive_orbits:
        candidates, sink_normalized_configuration_count = exhaustive_orbits()
    else:
        candidates = None
        sink_normalized_configuration_count = None
    target = len(candidates) if candidates is not None else args.samples
    seen: set[tuple[int, ...]] = set()
    best: list[dict[str, Any]] = []
    positive = 0
    while len(seen) < target:
        masks = (
            candidates[len(seen)]
            if candidates is not None
            else draw_candidate(rng)
        )
        if masks in seen:
            continue
        seen.add(masks)
        topology = complete_bipartite_rectangular_topology(
            masks, (1,) * 6, 2, 3, N
        )
        optimum = optimize_topology(
            topology,
            N,
            rng,
            args.starts,
            args.iterations,
            args.variable_grand,
        )
        if optimum is None:
            continue
        row = {
            "masks": list(masks),
            "divergence": [list(vector) for vector in topology.divergence],
            **optimum,
        }
        if float(row["objective"]) > 1e-8:
            positive += 1
        best.append(row)
        best.sort(key=lambda item: float(item["objective"]), reverse=True)
        del best[10:]
        if len(seen) % 10 == 0 or positive:
            print(
                json.dumps(
                    {
                        "tested": len(seen),
                        "positive": positive,
                        "best_objective": best[0]["objective"],
                    }
                ),
                flush=True,
            )
        if positive:
            break
    payload = {
        "status": (
            "positive_terminal_obstruction_found"
            if positive
            else "no_positive_terminal_obstruction_found"
        ),
        "n": N,
        "topology": "complete_bipartite_k23_all_source_pairs_and_sinks_sharp",
        "exhaustive_up_to_symmetry": args.exhaustive_orbits,
        "sink_normalized_configuration_count": (
            sink_normalized_configuration_count
        ),
        "orbit_count": None if candidates is None else len(candidates),
        "orbit_representatives": (
            None
            if candidates is None
            else [
                {
                    "masks": list(masks),
                    "source_disagreement_coordinates": (
                        source_disagreement_coordinates(masks)
                    ),
                }
                for masks in candidates
            ]
        ),
        "tested": len(seen),
        "positive": positive,
        "variable_grand": args.variable_grand,
        "configuration": vars(args) | {"output": str(args.output)},
        "best": best,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 1 if positive else 0


if __name__ == "__main__":
    raise SystemExit(main())
