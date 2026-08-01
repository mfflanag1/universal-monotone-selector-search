#!/usr/bin/env python3
"""Screen exact-family copies aligned by permutations and modular shifts."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import permutations
from pathlib import Path

from family_search import (
    common_core_gap_float,
    family_slack_float,
    serial_family,
)
from n5_intrinsic_permutation_union import (
    induced_family,
    permute_game,
)


F = Fraction


def affine_modular_transform(
    target: tuple[F, ...],
    source: tuple[F, ...],
    n: int,
    allow_scaling: bool,
) -> tuple[F, tuple[F, ...]] | None:
    def centered(game: tuple[F, ...], coalition: int) -> F:
        return game[coalition] - sum(
            (
                game[1 << player]
                for player in range(n)
                if coalition >> player & 1
            ),
            F(0),
        )

    scale = F(1)
    if allow_scaling:
        scale = F(0)
        for coalition in range(1 << n):
            source_value = centered(source, coalition)
            target_value = centered(target, coalition)
            if source_value == 0:
                if target_value != 0:
                    return None
                continue
            candidate = target_value / source_value
            if scale == 0:
                scale = candidate
            elif scale != candidate:
                return None
        if scale <= 0:
            return None
    coefficients = tuple(
        target[1 << player] - scale * source[1 << player]
        for player in range(n)
    )
    shift = tuple(
        sum(
            (
                coefficients[player]
                for player in range(n)
                if coalition >> player & 1
            ),
            F(0),
        )
        for coalition in range(1 << n)
    )
    if any(
        target[coalition] - scale * source[coalition]
        != shift[coalition]
        for coalition in range(1 << n)
    ):
        return None
    return scale, shift


def monotone(game: tuple[F, ...], n: int) -> bool:
    return all(
        game[coalition] <= game[coalition | (1 << player)]
        for coalition in range(1 << n)
        for player in range(n)
        if not coalition >> player & 1
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--allow-scaling", action="store_true")
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    n = int(payload["n"])
    base_games = [
        tuple(F(value) for value in game)
        for game in payload["family"]["games"]
    ]
    base_image = tuple(sorted(base_games))
    images: dict[
        tuple[tuple[F, ...], ...],
        dict[str, object],
    ] = {}
    for permutation in permutations(range(n)):
        permuted = [
            permute_game(game, permutation)
            for game in base_games
        ]
        for target_node, target in enumerate(base_games):
            for source_node, source in enumerate(permuted):
                transform = affine_modular_transform(
                    target, source, n, args.allow_scaling
                )
                if transform is None:
                    continue
                scale, shift = transform
                translated = tuple(
                    sorted(
                        tuple(
                            scale * value + shift[coalition]
                            for coalition, value in enumerate(game)
                        )
                        for game in permuted
                    )
                )
                if translated == base_image or any(
                    not monotone(game, n) for game in translated
                ):
                    continue
                images.setdefault(
                    translated,
                    {
                        "permutation": list(permutation),
                        "target_node": target_node,
                        "source_node": source_node,
                        "scale": str(scale),
                        "shift": [str(value) for value in shift],
                    },
                )

    base_family = induced_family(
        base_games, n, "n5_modular_alignment_base"
    )
    candidates = []
    for image, metadata in images.items():
        family = induced_family(
            base_games + list(image),
            n,
            "n5_modular_alignment_pair",
        )
        shared_nodes = 2 * len(base_games) - len(family.games)
        cross_edges = (
            len(family.edges)
            - 2 * len(base_family.edges)
        )
        candidates.append(
            (
                -cross_edges,
                len(family.games),
                image,
                metadata,
                family,
                shared_nodes,
                cross_edges,
            )
        )
    candidates.sort(key=lambda row: row[:2])
    if args.limit:
        candidates = candidates[: args.limit]

    rows = []
    best = None
    for index, (
        _,
        _,
        image,
        metadata,
        family,
        shared_nodes,
        cross_edges,
    ) in enumerate(candidates, start=1):
        margin, _ = family_slack_float(
            family.games, family.edges, n
        )
        common_gap = common_core_gap_float(family.games, n)
        row = {
            **metadata,
            "shared_nodes": shared_nodes,
            "cross_edges": cross_edges,
            "node_count": len(family.games),
            "edge_count": len(family.edges),
            "common_core_budget_gap_float": common_gap,
            "max_min_monotonicity_margin_float": margin,
        }
        rows.append(row)
        if best is None or margin < best[0]:
            best = (margin, image, metadata, family, common_gap, row)
            print(
                json.dumps(
                    {
                        "screened": index,
                        "candidate_count": len(candidates),
                        **row,
                    }
                ),
                flush=True,
            )

    if best is None:
        raise RuntimeError("no nontrivial monotone modular alignment")
    margin, _, metadata, family, common_gap, best_row = best
    result = {
        "status": "modular_alignment_screen_complete",
        "n": n,
        "source": str(args.input),
        "distinct_alignments": len(images),
        "screened_alignments": len(candidates),
        "best_alignment": metadata,
        "node_count": len(family.games),
        "edge_count": len(family.edges),
        "common_core_budget_gap_float": common_gap,
        "max_min_monotonicity_margin_float": margin,
        "best_row": best_row,
        "rows": rows,
        "family": serial_family(family, margin, common_gap),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "distinct_alignments",
                    "screened_alignments",
                    "best_alignment",
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
