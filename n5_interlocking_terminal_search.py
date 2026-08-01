#!/usr/bin/env python3
"""Search minimal interlocking signed-indicator terminal topologies."""

from __future__ import annotations

import argparse
import itertools
import json
import random
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from n5_facet_search import load_facets
from n5_terminal_topology_relaxation import reconstruct_exact_dual


F = Fraction


def members(coalition: int, n: int) -> list[int]:
    return [player for player in range(n) if coalition >> player & 1]


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
            index for index, coalition in enumerate(sources)
            if coalition >> player & 1
        ]
        sink_nodes = [
            index for index, coalition in enumerate(sinks)
            if coalition >> player & 1
        ]
        if len(source_nodes) != len(sink_nodes):
            raise ValueError("source and sink incidence do not balance")
        if len(source_nodes) == 1:
            supports[source_nodes[0]][sink_nodes[0]] |= 1 << player
        elif len(source_nodes) == 2:
            crossing = crossing_mask >> doubled_index & 1
            doubled_index += 1
            for source_node in range(2):
                sink_node = source_node ^ crossing
                supports[source_node][sink_node] |= 1 << player
    return (tuple(supports[0]), tuple(supports[1]))


def scalar_block_transport_exists(
    sources: tuple[int, ...], sinks: tuple[int, ...], n: int
) -> bool:
    # Variables are w_rq.  The exact equations have only
    # coefficients zero and one, so a floating feasibility solve is decisive
    # at the coarse filtering stage; every retained topology is later exactified.
    source_count = len(sources)
    sink_count = len(sinks)
    variable_count = source_count * sink_count
    rows: list[list[float]] = []
    rhs: list[float] = []
    for source in range(source_count):
        for player in range(n):
            row = [0.0] * variable_count
            for sink in range(sink_count):
                if sinks[sink] >> player & 1:
                    row[sink_count * source + sink] = 1.0
            rows.append(row)
            rhs.append(float(sources[source] >> player & 1))
    for sink in range(sink_count):
        row = [0.0] * variable_count
        for source in range(source_count):
            row[sink_count * source + sink] = 1.0
        rows.append(row)
        rhs.append(1.0)
    result = linprog(
        np.zeros(variable_count),
        A_eq=np.asarray(rows),
        b_eq=np.asarray(rhs),
        bounds=[(0.0, None)] * variable_count,
        method="highs",
    )
    return bool(result.success)


class TerminalLP:
    def __init__(self, n: int, source_count: int) -> None:
        self.n = n
        self.grand = (1 << n) - 1
        self.coalition_count = self.grand + 1
        self.source_count = source_count
        self.node_count = 2 * source_count
        self.variable_count = self.node_count * self.coalition_count
        if n == 5:
            self.facets, _ = load_facets()
        elif n == 6:
            from n6_mixed_terminal_search import load_facets as load_n6_facets

            self.facets = load_n6_facets()
        else:
            raise ValueError("complete exact-cone facets are available only for n=5 or n=6")
        self.base_rows: list[dict[int, float]] = []
        self.base_rhs: list[float] = []
        self.base_metadata: list[tuple[Any, ...]] = []
        self.equalities: list[dict[int, float]] = []
        self.equality_rhs: list[float] = []
        self.equality_metadata: list[tuple[Any, ...]] = []
        for node in range(self.node_count):
            self.equalities.append({self.column(node, 0): 1.0})
            self.equality_rhs.append(0.0)
            self.equality_metadata.append(("empty", node))
            self.equalities.append({self.column(node, self.grand): 1.0})
            self.equality_rhs.append(1.0)
            self.equality_metadata.append(("grand", node))
            for facet_index, facet in enumerate(self.facets):
                self.base_rows.append(
                    {
                        self.column(node, coalition): -float(coefficient)
                        for coalition, coefficient in enumerate(facet)
                        if coefficient
                    }
                )
                self.base_rhs.append(0.0)
                self.base_metadata.append(
                    ("exact_cone", node, facet_index)
                )
            for coalition in range(self.coalition_count):
                for player in range(n):
                    if coalition >> player & 1:
                        continue
                    successor = coalition | (1 << player)
                    self.base_rows.append(
                        {
                            self.column(node, coalition): 1.0,
                            self.column(node, successor): -1.0,
                        }
                    )
                    self.base_rhs.append(0.0)
                    self.base_metadata.append(
                        (
                            "game_monotonicity",
                            node,
                            coalition,
                            successor,
                        )
                    )
        self.base_equalities = self.matrix(
            self.equalities, self.variable_count
        )

    def column(self, node: int, coalition: int) -> int:
        return node * self.coalition_count + coalition

    @staticmethod
    def matrix(rows: list[dict[int, float]], columns: int):
        row_indices: list[int] = []
        column_indices: list[int] = []
        values: list[float] = []
        for row_index, row in enumerate(rows):
            for column, value in row.items():
                if value:
                    row_indices.append(row_index)
                    column_indices.append(column)
                    values.append(value)
        return coo_matrix(
            (values, (row_indices, column_indices)),
            shape=(len(rows), columns),
        ).tocsr()

    def solve(
        self,
        sources: tuple[int, ...],
        sinks: tuple[int, ...],
        supports: tuple[tuple[int, ...], ...],
        exactify: bool = False,
    ) -> tuple[float, np.ndarray, dict[str, Any] | None] | None:
        inequalities = list(self.base_rows)
        rhs = list(self.base_rhs)
        inequality_metadata = list(self.base_metadata)
        equalities = list(self.equalities)
        equality_rhs = list(self.equality_rhs)
        equality_metadata = list(self.equality_metadata)
        for source in range(self.source_count):
            for sink in range(self.source_count):
                protected = supports[source][sink]
                if not protected:
                    continue
                lower = source
                upper = self.source_count + sink
                for coalition in range(self.coalition_count):
                    row = {
                        self.column(lower, coalition): 1.0,
                        self.column(upper, coalition): -1.0,
                    }
                    if coalition & protected != protected:
                        equalities.append(row)
                        equality_rhs.append(0.0)
                        equality_metadata.append(
                            (
                                "protected_invariance",
                                source,
                                sink,
                                coalition,
                                protected,
                            )
                        )
                    else:
                        inequalities.append(row)
                        rhs.append(0.0)
                        inequality_metadata.append(
                            (
                                "directed_monotonicity",
                                source,
                                sink,
                                coalition,
                                protected,
                            )
                        )
        objective = np.zeros(self.variable_count)
        objective_exact = [F(0)] * self.variable_count
        for source, coalition in enumerate(sources):
            objective[self.column(source, coalition)] = -1.0
            objective_exact[self.column(source, coalition)] = F(-1)
        for sink, coalition in enumerate(sinks):
            column = self.column(
                self.source_count + sink, self.grand ^ coalition
            )
            objective[column] = -1.0
            objective_exact[column] = F(-1)
        result = linprog(
            objective,
            A_ub=self.matrix(inequalities, self.variable_count),
            b_ub=np.asarray(rhs),
            A_eq=self.matrix(equalities, self.variable_count),
            b_eq=np.asarray(equality_rhs),
            bounds=[(None, None)] * self.variable_count,
            method="highs",
        )
        if not result.success:
            return None
        # Each negative terminal contributes -G=-1 to lower expectation.
        exact_certificate = (
            reconstruct_exact_dual(
                result,
                inequalities,
                rhs,
                inequality_metadata,
                equalities,
                equality_rhs,
                equality_metadata,
                objective_exact,
                F(-self.source_count),
            )
            if exactify
            else None
        )
        return (
            -float(result.fun) - float(self.source_count),
            result.x,
            exact_certificate,
        )


def connected_incidence(supports: tuple[tuple[int, ...], ...]) -> bool:
    k = len(supports)
    adjacency = [set() for _ in range(2 * k)]
    for source in range(k):
        for sink in range(k):
            if supports[source][sink]:
                adjacency[source].add(k + sink)
                adjacency[k + sink].add(source)
    if any(len(neighbors) < 2 for neighbors in adjacency):
        return False
    seen = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for neighbor in adjacency[node]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return len(seen) == 2 * k


def random_routing(
    k: int, n: int, rng: random.Random
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[tuple[int, ...], ...]]:
    sources = [0] * k
    sinks = [0] * k
    supports = [[0] * k for _ in range(k)]
    for player in range(n):
        degree = rng.randint(1, k)
        source_nodes = rng.sample(range(k), degree)
        sink_nodes = rng.sample(range(k), degree)
        rng.shuffle(sink_nodes)
        for source, sink in zip(source_nodes, sink_nodes, strict=True):
            sources[source] |= 1 << player
            sinks[sink] |= 1 << player
            supports[source][sink] |= 1 << player
    return tuple(sources), tuple(sinks), tuple(
        tuple(row) for row in supports
    )


def canonical_two_by_two_routings(n: int):
    """One representative per player/source/sink orbit of interlocking routings."""

    def compositions(total: int, parts: int, prefix: tuple[int, ...] = ()):
        if parts == 1:
            yield prefix + (total,)
            return
        for value in range(total + 1):
            yield from compositions(total - value, parts - 1, prefix + (value,))

    def transform(counts: tuple[int, ...], source_swap: int, sink_swap: int):
        cells = [0] * 4
        for source in range(2):
            for sink in range(2):
                cells[2 * (source ^ source_swap) + (sink ^ sink_swap)] = counts[2 * source + sink]
        straight, crossing = counts[4:6]
        if source_swap ^ sink_swap:
            straight, crossing = crossing, straight
        return tuple(cells) + (straight, crossing, counts[6])

    representatives = {
        min(transform(counts, source_swap, sink_swap) for source_swap in range(2) for sink_swap in range(2))
        for counts in compositions(n, 7)
    }
    for counts in sorted(representatives):
        players = iter(range(n))
        sources = [0, 0]
        sinks = [0, 0]
        supports = [[0, 0], [0, 0]]
        for cell, count in enumerate(counts[:4]):
            source, sink = divmod(cell, 2)
            for _ in range(count):
                player = next(players)
                bit = 1 << player
                sources[source] |= bit
                sinks[sink] |= bit
                supports[source][sink] |= bit
        for crossing, count in enumerate(counts[4:6]):
            for _ in range(count):
                player = next(players)
                bit = 1 << player
                sources[0] |= bit
                sources[1] |= bit
                sinks[0] |= bit
                sinks[1] |= bit
                supports[0][crossing] |= bit
                supports[1][crossing ^ 1] |= bit
        for _ in range(counts[6]):
            next(players)
        if all(sources) and all(sinks) and all(
            supports[source][sink]
            for source in range(2)
            for sink in range(2)
        ):
            yield tuple(sources), tuple(sinks), None, tuple(tuple(row) for row in supports)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=5)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--samples", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--exactify", action="store_true")
    parser.add_argument("--canonical-types", action="store_true")
    parser.add_argument("--include-grand", action="store_true")
    args = parser.parse_args()
    n = args.n
    grand = (1 << n) - 1
    terminal_lp = TerminalLP(n, args.k)
    tested = 0
    outside_block_transport = 0
    positive = 0
    exactified = 0
    exactification_failures = 0
    exact_certificates = []
    best: list[dict[str, Any]] = []
    candidates: Any
    if args.k == 2:
        if args.canonical_types:
            candidates = canonical_two_by_two_routings(n)
        else:
            pairs_by_incidence: dict[
                tuple[int, ...], list[tuple[int, int]]
            ] = defaultdict(list)
            for left in range(1, grand):
                for right in range(left, grand):
                    pairs_by_incidence[incidence_sum((left, right), n)].append(
                        (left, right)
                    )
            exhaustive = []
            for pairs in pairs_by_incidence.values():
                for sources_index, sources in enumerate(pairs):
                    for sinks_index, sinks in enumerate(pairs):
                        if sinks_index == sources_index:
                            continue
                        doubled = sum(
                            incidence == 2
                            for incidence in incidence_sum(sources, n)
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
                                exhaustive.append(
                                    (sources, sinks, crossing_mask, supports)
                                )
            candidates = exhaustive
    else:
        rng = random.Random(args.seed)
        candidates = (
            (*random_routing(args.k, n, rng), None)
            for _ in range(args.samples)
        )

    for candidate in candidates:
        if args.k == 2:
            sources, sinks, crossing_mask, supports = candidate
        else:
            sources, sinks, supports, crossing_mask = candidate
        excluded = (0,) if args.include_grand else (0, grand)
        if any(coalition in excluded for coalition in (*sources, *sinks)):
            continue
        if args.k > 2 and not connected_incidence(supports):
            continue
        tested += 1
        if scalar_block_transport_exists(sources, sinks, n):
            continue
        outside_block_transport += 1
        solved = terminal_lp.solve(
            sources, sinks, supports, exactify=args.exactify
        )
        if solved is None:
            continue
        objective, point, exact_certificate = solved
        if objective > 1e-9:
            positive += 1
        if args.exactify:
            if exact_certificate is None:
                exactification_failures += 1
            else:
                exactified += 1
                exact_certificates.append(
                    {
                        "sources": list(sources),
                        "sinks": list(sinks),
                        "crossing_mask": crossing_mask,
                        "supports": [
                            list(support) for support in supports
                        ],
                        "exact_dual_certificate": exact_certificate,
                    }
                )
        row = {
            "sources": list(sources),
            "sinks": list(sinks),
            "crossing_mask": crossing_mask,
            "supports": [list(support) for support in supports],
            "maximum_lower_expectation_objective_float": objective,
            "terminal_games": [
                [
                    float(point[node * (grand + 1) + coalition])
                    for coalition in range(grand + 1)
                ]
                for node in range(2 * args.k)
            ],
            **(
                {}
                if exact_certificate is None
                else {"exact_dual_certificate": exact_certificate}
            ),
        }
        best.append(row)
        best.sort(
            key=lambda candidate: candidate[
                "maximum_lower_expectation_objective_float"
            ],
            reverse=True,
        )
        del best[args.top :]
        if args.limit is not None and tested >= args.limit:
            break

    result = {
        "status": (
            "interlocking_terminal_classification_exact"
            if args.exactify and exactification_failures == 0
            else "interlocking_terminal_search_float"
        ),
        "n": n,
        "source_and_sink_count": args.k,
        "tested_fully_interlocking_routings": tested,
        "outside_scalar_block_transport": outside_block_transport,
        "positive_terminal_relaxations": positive,
        "exactified_zero_bounds": exactified,
        "exactification_failures": exactification_failures,
        **(
            {}
            if not args.exactify
            else {"exact_certificates": exact_certificates}
        ),
        "best": best,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
