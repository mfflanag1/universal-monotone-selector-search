#!/usr/bin/env python3
"""Certify the five-player 2-1-2 nested lower-expectation identity."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
sys.path.insert(0, str(PROJECT_SRC))

from n5_facet_search import load_facets

from n5_mixed_divergence_premium import choquet_representation
from n5_mixed_terminal_branch_search import extreme_representations


F = Fraction


def coefficient_vector(
    representation: tuple[F, tuple[tuple[int, F], ...]], n: int
) -> tuple[F, ...]:
    grand = (1 << n) - 1
    beta, terms = representation
    result = [F(0)] * (grand + 1)
    result[grand] = beta
    for coalition, weight in terms:
        result[coalition] += weight
    return tuple(result)


def exact_generators(n: int) -> tuple[list[tuple[F, ...]], list[tuple[Any, ...]]]:
    grand = (1 << n) - 1
    facets, _ = load_facets()
    generators = [tuple(F(value) for value in facet) for facet in facets]
    metadata: list[tuple[Any, ...]] = [
        ("exact_facet", index) for index in range(len(facets))
    ]
    for sign in (-1, 1):
        empty_row = [F(0)] * (grand + 1)
        empty_row[0] = F(sign)
        generators.append(tuple(empty_row))
        metadata.append(("empty_equality", sign))
    for coalition in range(grand + 1):
        for player in range(n):
            if coalition >> player & 1:
                continue
            successor = coalition | (1 << player)
            row = [F(0)] * (grand + 1)
            row[coalition] = F(-1)
            row[successor] = F(1)
            generators.append(tuple(row))
            metadata.append(
                ("game_monotonicity", coalition, successor)
            )
    return generators, metadata


def find_certificate(
    target: tuple[F, ...], generators: list[tuple[F, ...]]
) -> list[tuple[int, F]]:
    if not any(target):
        return []
    matrix = csc_matrix(
        np.asarray(
            [
                [float(generator[row]) for generator in generators]
                for row in range(len(target))
            ]
        )
    )
    objective = np.asarray(
        [1.0 + (index % 97) * 1e-7 for index in range(len(generators))]
    )
    solved = linprog(
        objective,
        A_eq=matrix,
        b_eq=np.asarray([float(value) for value in target]),
        bounds=[(0.0, None)] * len(generators),
        method="highs-ds",
        options={
            "dual_feasibility_tolerance": 1e-10,
            "primal_feasibility_tolerance": 1e-10,
        },
    )
    if not solved.success:
        raise RuntimeError("target is outside the generated dual cone")
    support = [
        (index, F(float(weight)).limit_denominator(1_000_000))
        for index, weight in enumerate(solved.x)
        if weight > 1e-9
    ]
    reconstructed = [F(0)] * len(target)
    for index, weight in support:
        for coordinate, coefficient in enumerate(generators[index]):
            reconstructed[coordinate] += weight * coefficient
    if tuple(reconstructed) != target:
        raise RuntimeError("floating support did not rationalize exactly")
    if any(weight < 0 for _, weight in support):
        raise RuntimeError("certificate has a negative cone coefficient")
    return support


def generate(output: Path) -> None:
    n = 5
    divergence = (F(-1), F(-1), F(0), F(1), F(1))
    choquet = choquet_representation(divergence, n)
    choquet_vector = coefficient_vector(choquet, n)
    generators, metadata = exact_generators(n)
    certificates = []
    for branch, representation in enumerate(
        extreme_representations(list(divergence), n)
    ):
        branch_vector = coefficient_vector(representation, n)
        target = tuple(
            choquet_vector[index] - branch_vector[index]
            for index in range(1 << n)
        )
        support = find_certificate(target, generators)
        beta, terms = representation
        certificates.append(
            {
                "branch": branch,
                "beta": str(beta),
                "terms": [
                    [coalition, str(weight)] for coalition, weight in terms
                ],
                "dual_cone_support": [
                    {
                        "generator": list(metadata[index]),
                        "weight": str(weight),
                    }
                    for index, weight in support
                ],
            }
        )
    payload = {
        "status": "n5_middle_split_identity_exact_certificate",
        "n": n,
        "divergence": [str(value) for value in divergence],
        "choquet_beta": str(choquet[0]),
        "choquet_terms": [
            [coalition, str(weight)] for coalition, weight in choquet[1]
        ],
        "branch_count": len(certificates),
        "certificates": certificates,
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")


def verify(archive: Path) -> None:
    payload = json.loads(archive.read_text())
    if payload["status"] != "n5_middle_split_identity_exact_certificate":
        raise RuntimeError("unexpected archive status")
    n = int(payload["n"])
    divergence = tuple(F(value) for value in payload["divergence"])
    expected_representations = extreme_representations(list(divergence), n)
    if int(payload["branch_count"]) != len(expected_representations):
        raise RuntimeError("branch count mismatch")
    generators, metadata = exact_generators(n)
    generator_lookup = {
        tuple(row): index for index, row in enumerate(metadata)
    }
    choquet = (
        F(payload["choquet_beta"]),
        tuple(
            (int(coalition), F(weight))
            for coalition, weight in payload["choquet_terms"]
        ),
    )
    if choquet != choquet_representation(divergence, n):
        raise RuntimeError("Choquet representation mismatch")
    choquet_vector = coefficient_vector(choquet, n)
    seen = set()
    for entry in payload["certificates"]:
        branch = int(entry["branch"])
        if branch in seen or not 0 <= branch < len(expected_representations):
            raise RuntimeError("invalid or duplicate branch")
        seen.add(branch)
        representation = (
            F(entry["beta"]),
            tuple(
                (int(coalition), F(weight))
                for coalition, weight in entry["terms"]
            ),
        )
        if representation != expected_representations[branch]:
            raise RuntimeError("branch representation mismatch")
        branch_vector = coefficient_vector(representation, n)
        target = tuple(
            choquet_vector[index] - branch_vector[index]
            for index in range(1 << n)
        )
        reconstructed = [F(0)] * (1 << n)
        for active in entry["dual_cone_support"]:
            generator = tuple(active["generator"])
            if generator not in generator_lookup:
                raise RuntimeError("unknown dual-cone generator")
            weight = F(active["weight"])
            if weight < 0:
                raise RuntimeError("negative dual-cone weight")
            for coordinate, coefficient in enumerate(
                generators[generator_lookup[generator]]
            ):
                reconstructed[coordinate] += weight * coefficient
        if tuple(reconstructed) != target:
            raise RuntimeError("exact cone identity failed")
    if len(seen) != len(expected_representations):
        raise RuntimeError("archive omits a lower-expectation branch")
    print(
        "PASS: exact five-player 2-1-2 identity, "
        f"{len(seen)} lower-expectation branches"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--generate", action="store_true")
    args = parser.parse_args()
    if args.generate:
        generate(args.archive)
    verify(args.archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
