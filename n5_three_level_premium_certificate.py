#!/usr/bin/env python3
"""Certify sharp five-player three-level lower-expectation premiums."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from n5_middle_split_identity_certificate import (
    coefficient_vector,
    exact_generators,
    find_certificate,
)
from n5_mixed_divergence_premium import choquet_representation
from n5_mixed_terminal_branch_search import extreme_representations


F = Fraction


def dot(left: tuple[F, ...], right: tuple[F, ...]) -> F:
    return sum((a * b for a, b in zip(left, right, strict=True)), F(0))


def generate(screen: Path, output: Path, four_level: bool = False) -> None:
    n = 5
    grand = (1 << n) - 1
    screen_payload = json.loads(screen.read_text())
    generators, metadata = exact_generators(n)
    results: list[dict[str, Any]] = []
    for screened in screen_payload["results"]:
        divergence = tuple(F(value) for value in screened["divergence"])
        premium = F(float(screened["premium"])).limit_denominator(1_000_000)
        choquet = choquet_representation(divergence, n)
        choquet_vector = coefficient_vector(choquet, n)
        certificates = []
        representations = extreme_representations(list(divergence), n)
        for branch, representation in enumerate(representations):
            branch_vector = coefficient_vector(representation, n)
            target = tuple(
                choquet_vector[index]
                - branch_vector[index]
                + (premium if index == grand else F(0))
                for index in range(1 << n)
            )
            support = find_certificate(target, generators)
            beta, terms = representation
            certificates.append(
                {
                    "branch": branch,
                    "beta": str(beta),
                    "terms": [
                        [coalition, str(weight)]
                        for coalition, weight in terms
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
        witness_game = tuple(
            F(float(value)).limit_denominator(1_000_000)
            for value in screened["game"]
        )
        witness_representation = (
            representations[int(screened["branch"])]
            if "branch" in screened
            else (
                F(screened["beta"]),
                tuple(
                    (int(coalition), F(weight))
                    for coalition, weight in screened["terms"]
                ),
            )
        )
        if witness_representation not in representations:
            raise RuntimeError("screened witness branch is not extreme")
        results.append(
            {
                "level_multiplicities": screened["level_multiplicities"],
                "divergence": [str(value) for value in divergence],
                "premium": str(premium),
                "choquet_beta": str(choquet[0]),
                "choquet_terms": [
                    [coalition, str(weight)]
                    for coalition, weight in choquet[1]
                ],
                "branch_count": len(certificates),
                "certificates": certificates,
                "witness_game": [str(value) for value in witness_game],
                "witness_branch": representations.index(
                    witness_representation
                ),
            }
        )
    payload = {
        "status": (
            "n5_four_level_premium_exact_certificate"
            if four_level
            else "n5_three_level_premium_exact_certificate"
        ),
        "n": n,
        "results": results,
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")


def verify(archive: Path) -> None:
    payload = json.loads(archive.read_text())
    if payload["status"] not in {
        "n5_three_level_premium_exact_certificate",
        "n5_four_level_premium_exact_certificate",
    }:
        raise RuntimeError("unexpected archive status")
    n = int(payload["n"])
    grand = (1 << n) - 1
    generators, metadata = exact_generators(n)
    generator_lookup = {
        tuple(row): index for index, row in enumerate(metadata)
    }
    expected_multiplicities = (
        {(1, 1, 1, 2), (1, 1, 2, 1), (1, 2, 1, 1), (2, 1, 1, 1)}
        if payload["status"] == "n5_four_level_premium_exact_certificate"
        else {
            (negative, middle, n - negative - middle)
            for negative in range(1, n - 1)
            for middle in range(1, n - negative)
        }
    )
    seen_multiplicities = set()
    total_branches = 0
    for result in payload["results"]:
        multiplicities = tuple(int(value) for value in result["level_multiplicities"])
        if multiplicities in seen_multiplicities:
            raise RuntimeError("duplicate level multiplicities")
        seen_multiplicities.add(multiplicities)
        divergence = tuple(F(value) for value in result["divergence"])
        premium = F(result["premium"])
        choquet = (
            F(result["choquet_beta"]),
            tuple(
                (int(coalition), F(weight))
                for coalition, weight in result["choquet_terms"]
            ),
        )
        if choquet != choquet_representation(divergence, n):
            raise RuntimeError("Choquet representation mismatch")
        choquet_vector = coefficient_vector(choquet, n)
        representations = extreme_representations(list(divergence), n)
        if int(result["branch_count"]) != len(representations):
            raise RuntimeError("branch count mismatch")
        seen_branches = set()
        for entry in result["certificates"]:
            branch = int(entry["branch"])
            if branch in seen_branches or not 0 <= branch < len(representations):
                raise RuntimeError("invalid or duplicate branch")
            seen_branches.add(branch)
            representation = (
                F(entry["beta"]),
                tuple(
                    (int(coalition), F(weight))
                    for coalition, weight in entry["terms"]
                ),
            )
            if representation != representations[branch]:
                raise RuntimeError("branch representation mismatch")
            branch_vector = coefficient_vector(representation, n)
            target = tuple(
                choquet_vector[index]
                - branch_vector[index]
                + (premium if index == grand else F(0))
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
                raise RuntimeError("exact upper-bound identity failed")
        if len(seen_branches) != len(representations):
            raise RuntimeError("archive omits a lower-expectation branch")

        witness = tuple(F(value) for value in result["witness_game"])
        if witness[0] != 0 or witness[grand] != 1:
            raise RuntimeError("witness normalization failed")
        if any(dot(generator, witness) < 0 for generator in generators):
            raise RuntimeError("witness is not a monotone exact game")
        witness_representation = representations[int(result["witness_branch"])]
        witness_gap = dot(
            tuple(
                a - b
                for a, b in zip(
                    coefficient_vector(witness_representation, n),
                    choquet_vector,
                    strict=True,
                )
            ),
            witness,
        )
        if witness_gap != premium:
            raise RuntimeError("witness does not attain the sharp premium")
        total_branches += len(representations)
    if seen_multiplicities != expected_multiplicities:
        raise RuntimeError("archive omits a three-level multiplicity type")
    print(
        "PASS: exact sharp five-player mixed-divergence premiums, "
        f"{total_branches} lower-expectation branches"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--generate-from", type=Path)
    parser.add_argument("--four-level", action="store_true")
    args = parser.parse_args()
    if args.generate_from:
        generate(args.generate_from, args.archive, args.four_level)
    verify(args.archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
