#!/usr/bin/env python3
"""Audit protection locality of a paired-slack exact-cover lift."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np


F = Fraction


def core_vertices(game: list[F], n: int) -> list[tuple[F, ...]]:
    grand = (1 << n) - 1
    vertices: set[tuple[F, ...]] = set()
    for tight in itertools.combinations(range(1, grand), n - 1):
        matrix = [[F(1) for _ in range(n)]]
        rhs = [game[grand]]
        for coalition in tight:
            matrix.append(
                [F(bool(coalition >> player & 1)) for player in range(n)]
            )
            rhs.append(game[coalition])
        augmented = [row + [value] for row, value in zip(matrix, rhs, strict=True)]
        rank = 0
        for column in range(n):
            pivot = next(
                (row for row in range(rank, n) if augmented[row][column]),
                None,
            )
            if pivot is None:
                continue
            augmented[rank], augmented[pivot] = augmented[pivot], augmented[rank]
            scale = augmented[rank][column]
            augmented[rank] = [value / scale for value in augmented[rank]]
            for row in range(n):
                if row == rank or not augmented[row][column]:
                    continue
                scale = augmented[row][column]
                augmented[row] = [
                    left - scale * right
                    for left, right in zip(
                        augmented[row], augmented[rank], strict=True
                    )
                ]
            rank += 1
        if rank != n:
            continue
        point = tuple(augmented[player][-1] for player in range(n))
        if all(
            sum(
                point[player]
                for player in range(n)
                if coalition >> player & 1
            )
            >= game[coalition]
            for coalition in range(1, grand)
        ):
            vertices.add(point)
    if not vertices:
        raise RuntimeError("core has no enumerated vertices")
    return sorted(vertices)


def lower(vertices: list[tuple[F, ...]], direction: tuple[int, ...]) -> F:
    return min(
        sum(F(coefficient) * value for coefficient, value in zip(direction, point, strict=True))
        for point in vertices
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = json.loads(args.source.read_text())["family"]
    games = [[F(value) for value in game] for game in raw["games"]]
    n = (len(games[0]) - 1).bit_length()
    grand = (1 << n) - 1
    constraints = tuple(range(1, grand))
    vertices = [core_vertices(game, n) for game in games]
    cache: list[dict[tuple[int, ...], F]] = [dict() for _ in games]

    def lifted_worth(game_index: int, original: int, negative_helpers: int) -> F:
        direction = [int(bool(original >> player & 1)) for player in range(n)]
        constant = F(0)
        for index, coalition in enumerate(constraints):
            if not negative_helpers >> index & 1:
                continue
            constant += games[game_index][coalition]
            for player in range(n):
                if coalition >> player & 1:
                    direction[player] -= 1
        key = tuple(direction)
        value = cache[game_index].get(key)
        if value is None:
            value = lower(vertices[game_index], key)
            cache[game_index][key] = value
        return constant + value

    rows = []
    first_violation = None
    for edge_index, edge in enumerate(raw["edges"]):
        lower_game = int(edge["lower"])
        upper_game = int(edge["upper"])
        protected = int(edge["coalition"])
        changed = 0
        leaking = 0
        minimum_delta = None
        maximum_delta = F(0)
        for original in range(grand + 1):
            for negative_helpers in range(1 << len(constraints)):
                low = lifted_worth(lower_game, original, negative_helpers)
                high = lifted_worth(upper_game, original, negative_helpers)
                delta = high - low
                minimum_delta = delta if minimum_delta is None else min(minimum_delta, delta)
                maximum_delta = max(maximum_delta, delta)
                if delta <= 0:
                    continue
                changed += 1
                if original & protected != protected:
                    leaking += 1
                    if first_violation is None:
                        first_violation = {
                            "edge": edge_index,
                            "protected": protected,
                            "original_members": original,
                            "negative_helper_constraints": [
                                coalition
                                for index, coalition in enumerate(constraints)
                                if negative_helpers >> index & 1
                            ],
                            "delta": str(delta),
                        }
        rows.append(
            {
                "edge": edge_index,
                "protected": protected,
                "changed_patterns": changed,
                "leaking_patterns": leaking,
                "minimum_delta": str(minimum_delta),
                "maximum_delta": str(maximum_delta),
            }
        )
        print(json.dumps(rows[-1]), flush=True)
    payload = {
        "status": "protection_local" if first_violation is None else "protection_leak",
        "n": n,
        "constraint_count": len(constraints),
        "core_vertex_counts": [len(points) for points in vertices],
        "distinct_direction_counts": [len(values) for values in cache],
        "rows": rows,
        "first_violation": first_violation,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0 if first_violation is None else 1


if __name__ == "__main__":
    raise SystemExit(main())
