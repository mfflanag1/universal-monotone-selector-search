#!/usr/bin/env python3
"""Independently verify the n=6 two-source/two-sink exact classification."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from sympy import Matrix

from n6_mixed_terminal_search import load_facets


F = Fraction
N = 6
GRAND = 63


def incidence(pair: tuple[int, int]) -> tuple[int, ...]:
    return tuple(sum(mask >> player & 1 for mask in pair) for player in range(N))


def route_supports(
    sources: tuple[int, int], sinks: tuple[int, int], crossing_mask: int
) -> tuple[tuple[int, int], tuple[int, int]]:
    supports = [[0, 0], [0, 0]]
    doubled = 0
    for player in range(N):
        source_nodes = [index for index, mask in enumerate(sources) if mask >> player & 1]
        sink_nodes = [index for index, mask in enumerate(sinks) if mask >> player & 1]
        if len(source_nodes) != len(sink_nodes):
            raise RuntimeError("unbalanced incidence")
        if len(source_nodes) == 1:
            supports[source_nodes[0]][sink_nodes[0]] |= 1 << player
        elif len(source_nodes) == 2:
            crossing = crossing_mask >> doubled & 1
            doubled += 1
            supports[0][crossing] |= 1 << player
            supports[1][crossing ^ 1] |= 1 << player
    return tuple(tuple(row) for row in supports)


def signature(
    sources: tuple[int, int],
    sinks: tuple[int, int],
    supports: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[int, ...]:
    counts = [0] * 7
    for player in range(N):
        bit = 1 << player
        cells = [
            2 * source + sink
            for source in range(2)
            for sink in range(2)
            if supports[source][sink] & bit
        ]
        if not cells:
            counts[6] += 1
        elif len(cells) == 1:
            counts[cells[0]] += 1
        elif cells == [0, 3]:
            counts[4] += 1
        elif cells == [1, 2]:
            counts[5] += 1
        else:
            raise RuntimeError("supports do not route each player bijectively")
        if sum(bool(mask & bit) for mask in sources) != len(cells):
            raise RuntimeError("source incidence disagrees with supports")
        if sum(bool(mask & bit) for mask in sinks) != len(cells):
            raise RuntimeError("sink incidence disagrees with supports")

    def transformed(source_swap: int, sink_swap: int) -> tuple[int, ...]:
        cells = [0] * 4
        for source in range(2):
            for sink in range(2):
                cells[2 * (source ^ source_swap) + (sink ^ sink_swap)] = counts[2 * source + sink]
        straight, crossing = counts[4:6]
        if source_swap ^ sink_swap:
            straight, crossing = crossing, straight
        return tuple(cells) + (straight, crossing, counts[6])

    return min(transformed(source_swap, sink_swap) for source_swap in range(2) for sink_swap in range(2))


def enumerate_signatures() -> dict[tuple[int, ...], tuple[tuple[int, int], tuple[int, int]]]:
    pairs_by_incidence: dict[tuple[int, ...], list[tuple[int, int]]] = defaultdict(list)
    for left in range(1, GRAND + 1):
        for right in range(left, GRAND + 1):
            pair = (left, right)
            pairs_by_incidence[incidence(pair)].append(pair)
    representatives = {}
    for pairs in pairs_by_incidence.values():
        for sources in pairs:
            for sinks in pairs:
                doubled = sum(value == 2 for value in incidence(sources))
                for crossing_mask in range(1 << doubled):
                    supports = route_supports(sources, sinks, crossing_mask)
                    if all(supports[source][sink] for source in range(2) for sink in range(2)):
                        representatives.setdefault(signature(sources, sinks, supports), (sources, sinks))
    return representatives


def scalar_transport_exists(sources: tuple[int, int], sinks: tuple[int, int]) -> bool:
    rows = []
    rhs = []
    for source in range(2):
        for player in range(N):
            rows.append(
                [
                    int(column // 2 == source) * int(sinks[column % 2] >> player & 1)
                    for column in range(4)
                ]
            )
            rhs.append(int(sources[source] >> player & 1))
    for sink in range(2):
        rows.append([int(column % 2 == sink) for column in range(4)])
        rhs.append(1)
    matrix = Matrix(rows)
    target = Matrix(rhs)
    rank = matrix.rank()
    for columns in itertools.combinations(range(4), rank):
        basis = matrix[:, columns]
        if basis.rank() != rank or basis.row_join(target).rank() != rank:
            continue
        solution = list(basis.gauss_jordan_solve(target)[0])
        if all(value >= 0 for value in solution):
            return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("classification", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.classification.read_text())
    if payload["status"] != "interlocking_terminal_classification_exact":
        raise RuntimeError("classification is not exact")
    expected = enumerate_signatures()
    if len(expected) != 103:
        raise RuntimeError(f"independent orbit count is {len(expected)}, not 103")
    nontransport = {
        key for key, (sources, sinks) in expected.items()
        if not scalar_transport_exists(sources, sinks)
    }
    if len(nontransport) != 50:
        raise RuntimeError(f"independent nontransport count is {len(nontransport)}, not 50")

    facets = load_facets()
    coalition_count = GRAND + 1
    variable_count = 4 * coalition_count

    def column(node: int, coalition: int) -> int:
        return node * coalition_count + coalition

    archived = set()
    inequality_count = 0
    equality_count = 0
    for entry in payload["exact_certificates"]:
        sources = tuple(int(value) for value in entry["sources"])
        sinks = tuple(int(value) for value in entry["sinks"])
        supports = tuple(tuple(int(value) for value in row) for row in entry["supports"])
        key = signature(sources, sinks, supports)
        if key not in nontransport or key in archived:
            raise RuntimeError("unexpected or duplicate certificate topology")
        archived.add(key)
        objective = [F(0)] * variable_count
        objective[column(0, sources[0])] -= 1
        objective[column(1, sources[1])] -= 1
        objective[column(2, GRAND ^ sinks[0])] -= 1
        objective[column(3, GRAND ^ sinks[1])] -= 1
        stationarity = [F(0)] * variable_count
        dual_objective = F(0)

        def add(row: dict[int, F], weight: F) -> None:
            for coordinate, coefficient in row.items():
                stationarity[coordinate] += weight * coefficient

        certificate = entry["exact_dual_certificate"]
        if F(certificate["variable_objective"]) != -2:
            raise RuntimeError("wrong variable objective")
        if F(certificate["lower_expectation_objective_with_constant"]) != 0:
            raise RuntimeError("terminal bound is not zero")
        seen = set()
        for active in certificate["active_inequalities"]:
            metadata = tuple(active["row"])
            if metadata in seen:
                raise RuntimeError("duplicate inequality row")
            seen.add(metadata)
            weight = F(active["weight"])
            if weight >= 0:
                raise RuntimeError("inequality multiplier has wrong sign")
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
                node, coalition, successor = map(int, metadata[1:])
                if coalition & successor != coalition or (successor ^ coalition).bit_count() != 1:
                    raise RuntimeError("invalid game-monotonicity row")
                add({column(node, coalition): F(1), column(node, successor): F(-1)}, weight)
            elif row_type == "directed_monotonicity":
                _, source, sink, coalition, protected = metadata
                source, sink, coalition, protected = map(int, (source, sink, coalition, protected))
                if protected != supports[source][sink] or coalition & protected != protected:
                    raise RuntimeError("invalid directed row")
                add({column(source, coalition): F(1), column(2 + sink, coalition): F(-1)}, weight)
            else:
                raise RuntimeError(f"unknown inequality row {row_type}")
            inequality_count += 1
        for active in certificate["active_equalities"]:
            metadata = tuple(active["row"])
            if metadata in seen:
                raise RuntimeError("duplicate equality row")
            seen.add(metadata)
            weight = F(active["weight"])
            row_type = metadata[0]
            if row_type == "empty":
                add({column(int(metadata[1]), 0): F(1)}, weight)
            elif row_type == "grand":
                add({column(int(metadata[1]), GRAND): F(1)}, weight)
                dual_objective += weight
            elif row_type == "protected_invariance":
                _, source, sink, coalition, protected = metadata
                source, sink, coalition, protected = map(int, (source, sink, coalition, protected))
                if protected != supports[source][sink] or coalition & protected == protected:
                    raise RuntimeError("invalid invariance row")
                add({column(source, coalition): F(1), column(2 + sink, coalition): F(-1)}, weight)
            else:
                raise RuntimeError(f"unknown equality row {row_type}")
            equality_count += 1
        if stationarity != objective or dual_objective != -2:
            raise RuntimeError("rational dual certificate failed")

    if archived != nontransport:
        raise RuntimeError(f"archive misses {len(nontransport - archived)} nontransport orbits")
    expected_counts = {
        "tested_fully_interlocking_routings": 103,
        "outside_scalar_block_transport": 50,
        "positive_terminal_relaxations": 0,
        "exactified_zero_bounds": 50,
        "exactification_failures": 0,
    }
    for name, expected_value in expected_counts.items():
        if int(payload[name]) != expected_value:
            raise RuntimeError(f"reported {name} is inconsistent")
    additive = [F(coalition.bit_count(), N) for coalition in range(coalition_count)]
    for facet in facets:
        if sum((F(coefficient) * additive[index] for index, coefficient in enumerate(facet)), F(0)) < 0:
            raise RuntimeError("common additive primal violates exactness")
    print(
        "PASS: 103 n=6 interlocking k=2 orbits; 53 exact block transports and "
        f"50 rational zero-bound duals ({inequality_count} inequality, {equality_count} equality rows)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
