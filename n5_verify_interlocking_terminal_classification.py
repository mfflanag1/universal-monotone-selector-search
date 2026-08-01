#!/usr/bin/env python3
"""Verify the exhaustive exact two-source/two-sink terminal classification."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

from sympy import Matrix

from n5_facet_search import load_facets


F = Fraction


def incidence_sum(pair: tuple[int, int], n: int) -> tuple[int, ...]:
    return tuple(
        int(pair[0] >> player & 1) + int(pair[1] >> player & 1)
        for player in range(n)
    )


def route_supports(
    sources: tuple[int, int],
    sinks: tuple[int, int],
    crossing_mask: int,
    n: int,
) -> tuple[tuple[int, int], tuple[int, int]]:
    supports = [[0, 0], [0, 0]]
    doubled_index = 0
    for player in range(n):
        source_nodes = [
            index
            for index, coalition in enumerate(sources)
            if coalition >> player & 1
        ]
        sink_nodes = [
            index
            for index, coalition in enumerate(sinks)
            if coalition >> player & 1
        ]
        if len(source_nodes) != len(sink_nodes):
            raise RuntimeError("unbalanced player incidence")
        if len(source_nodes) == 1:
            supports[source_nodes[0]][sink_nodes[0]] |= 1 << player
        elif len(source_nodes) == 2:
            crossing = crossing_mask >> doubled_index & 1
            doubled_index += 1
            for source_node in range(2):
                supports[source_node][source_node ^ crossing] |= 1 << player
    return tuple(tuple(row) for row in supports)


def exhaustive_topologies(
    n: int,
) -> set[
    tuple[
        tuple[int, int],
        tuple[int, int],
        int,
        tuple[tuple[int, int], tuple[int, int]],
    ]
]:
    grand = (1 << n) - 1
    pairs_by_incidence: dict[
        tuple[int, ...], list[tuple[int, int]]
    ] = defaultdict(list)
    for left in range(1, grand):
        for right in range(left, grand):
            pair = (left, right)
            pairs_by_incidence[incidence_sum(pair, n)].append(pair)
    topologies = set()
    for pairs in pairs_by_incidence.values():
        for source_index, sources in enumerate(pairs):
            for sink_index, sinks in enumerate(pairs):
                if sink_index == source_index:
                    continue
                doubled = sum(
                    value == 2 for value in incidence_sum(sources, n)
                )
                for crossing_mask in range(1 << doubled):
                    supports = route_supports(
                        sources, sinks, crossing_mask, n
                    )
                    if all(
                        supports[source][sink]
                        for source in range(2)
                        for sink in range(2)
                    ):
                        topologies.add(
                            (
                                sources,
                                sinks,
                                crossing_mask,
                                supports,
                            )
                        )
    return topologies


def scalar_transport_is_inconsistent(
    sources: tuple[int, int],
    sinks: tuple[int, int],
    n: int,
) -> bool:
    matrix = []
    rhs = []
    for source in range(2):
        for player in range(n):
            matrix.append(
                [
                    int(column // 2 == source)
                    * int(sinks[column % 2] >> player & 1)
                    for column in range(4)
                ]
            )
            rhs.append(int(sources[source] >> player & 1))
    for sink in range(2):
        matrix.append(
            [int(column % 2 == sink) for column in range(4)]
        )
        rhs.append(1)
    rank = Matrix(matrix).rank()
    augmented_rank = Matrix(
        [row + [value] for row, value in zip(matrix, rhs, strict=True)]
    ).rank()
    return augmented_rank > rank


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("classification", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.classification.read_text())
    if payload["status"] != "interlocking_terminal_classification_exact":
        raise RuntimeError("unexpected classification status")
    n = int(payload["n"])
    if n != 5 or int(payload["source_and_sink_count"]) != 2:
        raise RuntimeError("classification has unexpected dimensions")
    expected_topologies = exhaustive_topologies(n)
    if len(expected_topologies) != 750:
        raise RuntimeError("independent enumeration count changed")

    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_count = 4
    variable_count = node_count * coalition_count
    facets, _ = load_facets()

    def column(node: int, coalition: int) -> int:
        if not 0 <= node < node_count or not 0 <= coalition <= grand:
            raise RuntimeError("invalid game coordinate")
        return node * coalition_count + coalition

    archived_topologies = set()
    total_inequality_rows = 0
    total_equality_rows = 0
    for entry in payload["exact_certificates"]:
        sources = tuple(int(value) for value in entry["sources"])
        sinks = tuple(int(value) for value in entry["sinks"])
        crossing_mask = int(entry["crossing_mask"])
        supports = tuple(
            tuple(int(value) for value in row)
            for row in entry["supports"]
        )
        topology = (sources, sinks, crossing_mask, supports)
        if topology in archived_topologies:
            raise RuntimeError("duplicate archived topology")
        archived_topologies.add(topology)
        if topology not in expected_topologies:
            raise RuntimeError("archive contains a non-enumerated topology")
        if route_supports(sources, sinks, crossing_mask, n) != supports:
            raise RuntimeError("archived routing support is incorrect")
        if not scalar_transport_is_inconsistent(sources, sinks, n):
            raise RuntimeError("routing admits scalar block transport")

        objective = [F(0)] * variable_count
        objective[column(0, sources[0])] = F(-1)
        objective[column(1, sources[1])] = F(-1)
        objective[column(2, grand ^ sinks[0])] = F(-1)
        objective[column(3, grand ^ sinks[1])] = F(-1)
        stationarity = [F(0)] * variable_count
        dual_objective = F(0)
        certificate = entry["exact_dual_certificate"]
        if F(certificate["variable_objective"]) != -2:
            raise RuntimeError("reported variable objective is incorrect")
        if F(certificate["lower_expectation_objective_with_constant"]) != 0:
            raise RuntimeError("reported terminal objective is not zero")

        def add(row: dict[int, F], weight: F) -> None:
            for coordinate, coefficient in row.items():
                stationarity[coordinate] += weight * coefficient

        seen_rows = set()
        for active in certificate["active_inequalities"]:
            metadata = tuple(active["row"])
            if metadata in seen_rows:
                raise RuntimeError("duplicate inequality dual row")
            seen_rows.add(metadata)
            weight = F(active["weight"])
            if weight >= 0:
                raise RuntimeError("inequality dual weight is not negative")
            row_type = metadata[0]
            if row_type == "exact_cone":
                _, node, facet_index = metadata
                facet = facets[int(facet_index)]
                add(
                    {
                        column(int(node), coalition): -F(coefficient)
                        for coalition, coefficient in enumerate(facet)
                        if coefficient
                    },
                    weight,
                )
            elif row_type == "game_monotonicity":
                _, node, coalition, successor = metadata
                coalition = int(coalition)
                successor = int(successor)
                if coalition & successor != coalition or (
                    successor ^ coalition
                ).bit_count() != 1:
                    raise RuntimeError("invalid game-monotonicity row")
                add(
                    {
                        column(int(node), coalition): F(1),
                        column(int(node), successor): F(-1),
                    },
                    weight,
                )
            elif row_type == "directed_monotonicity":
                _, source, sink, coalition, protected = metadata
                source = int(source)
                sink = int(sink)
                coalition = int(coalition)
                protected = int(protected)
                if (
                    protected != supports[source][sink]
                    or coalition & protected != protected
                ):
                    raise RuntimeError("invalid directed-monotonicity row")
                add(
                    {
                        column(source, coalition): F(1),
                        column(2 + sink, coalition): F(-1),
                    },
                    weight,
                )
            else:
                raise RuntimeError(f"unsupported inequality row {row_type}")
            total_inequality_rows += 1

        for active in certificate["active_equalities"]:
            metadata = tuple(active["row"])
            if metadata in seen_rows:
                raise RuntimeError("duplicate equality dual row")
            seen_rows.add(metadata)
            weight = F(active["weight"])
            row_type = metadata[0]
            if row_type == "empty":
                _, node = metadata
                add({column(int(node), 0): F(1)}, weight)
            elif row_type == "grand":
                _, node = metadata
                add({column(int(node), grand): F(1)}, weight)
                dual_objective += weight
            elif row_type == "protected_invariance":
                _, source, sink, coalition, protected = metadata
                source = int(source)
                sink = int(sink)
                coalition = int(coalition)
                protected = int(protected)
                if (
                    protected != supports[source][sink]
                    or coalition & protected == protected
                ):
                    raise RuntimeError("invalid protected-invariance row")
                add(
                    {
                        column(source, coalition): F(1),
                        column(2 + sink, coalition): F(-1),
                    },
                    weight,
                )
            else:
                raise RuntimeError(f"unsupported equality row {row_type}")
            total_equality_rows += 1

        if stationarity != objective:
            raise RuntimeError("exact dual stationarity failed")
        if dual_objective != -2:
            raise RuntimeError("exact dual objective failed")

    if archived_topologies != expected_topologies:
        missing = expected_topologies - archived_topologies
        raise RuntimeError(f"classification omits {len(missing)} topologies")
    if int(payload["tested_fully_interlocking_routings"]) != 750:
        raise RuntimeError("reported topology count is inconsistent")
    if int(payload["outside_scalar_block_transport"]) != 750:
        raise RuntimeError("reported transport count is inconsistent")
    if int(payload["exactified_zero_bounds"]) != 750:
        raise RuntimeError("reported exactification count is inconsistent")
    if int(payload["exactification_failures"]) != 0:
        raise RuntimeError("classification reports exactification failures")
    if int(payload["positive_terminal_relaxations"]) != 0:
        raise RuntimeError("classification reports a positive relaxation")

    additive = [
        F(coalition.bit_count(), n) for coalition in range(coalition_count)
    ]
    if any(
        sum(
            (
                F(coefficient) * additive[coalition]
                for coalition, coefficient in enumerate(facet)
            ),
            F(0),
        )
        < 0
        for facet in facets
    ):
        raise RuntimeError("common additive primal violates exactness")
    print(
        "PASS: exhaustive 750-topology exact classification, "
        f"{total_inequality_rows} inequality dual rows, "
        f"{total_equality_rows} equality dual rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
