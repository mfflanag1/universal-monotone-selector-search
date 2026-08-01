#!/usr/bin/env python3
"""Test whether an indicator-flow dual has a block-transport decomposition."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from math import lcm
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from exact_lp import exact_linprog


F = Fraction


def grouped_algebraic_solve(
    sources: list[dict[str, Any]],
    sinks: list[dict[str, Any]],
    n: int,
    orientation: str,
) -> bool:
    source_mass: dict[int, F] = defaultdict(F)
    sink_mass: dict[int, F] = defaultdict(F)
    for source in sources:
        source_mass[int(source["coalition"])] += F(source["coefficient"])
    for sink in sinks:
        sink_mass[int(sink["coalition"])] += F(sink["coefficient"])
    source_types = sorted(source_mass)
    sink_types = sorted(sink_mass)
    variable_count = len(source_types) * len(sink_types)
    equalities: list[list[F]] = []
    rhs: list[F] = []
    if orientation == "fork":
        for source_index, source_type in enumerate(source_types):
            for player in range(n):
                row = [F(0)] * variable_count
                for sink_index, sink_type in enumerate(sink_types):
                    if sink_type >> player & 1:
                        row[source_index * len(sink_types) + sink_index] = F(1)
                equalities.append(row)
                rhs.append(
                    source_mass[source_type]
                    if source_type >> player & 1
                    else F(0)
                )
        for sink_index, sink_type in enumerate(sink_types):
            row = [F(0)] * variable_count
            for source_index in range(len(source_types)):
                row[source_index * len(sink_types) + sink_index] = F(1)
            equalities.append(row)
            rhs.append(sink_mass[sink_type])
    else:
        for source_index, source_type in enumerate(source_types):
            row = [F(0)] * variable_count
            for sink_index in range(len(sink_types)):
                row[source_index * len(sink_types) + sink_index] = F(1)
            equalities.append(row)
            rhs.append(source_mass[source_type])
        for sink_index, sink_type in enumerate(sink_types):
            for player in range(n):
                row = [F(0)] * variable_count
                for source_index, source_type in enumerate(source_types):
                    if source_type >> player & 1:
                        row[source_index * len(sink_types) + sink_index] = F(1)
                equalities.append(row)
                rhs.append(
                    sink_mass[sink_type]
                    if sink_type >> player & 1
                    else F(0)
                )
    nonnegative_rows = []
    for variable in range(variable_count):
        row = [F(0)] * variable_count
        row[variable] = F(-1)
        nonnegative_rows.append(row)
    status, _, _ = exact_linprog(
        [F(0)] * variable_count,
        nonnegative_rows,
        [F(0)] * variable_count,
        equalities,
        rhs,
    )
    return status == "optimal"


def transport_solve(
    variable_count: int,
    equalities: list[list[F]],
    rhs: list[F],
    forbidden: set[int],
) -> tuple[bool, list[F] | None, str]:
    row_indices = []
    columns = []
    values = []
    for row_index, row in enumerate(equalities):
        for column, value in enumerate(row):
            if value:
                row_indices.append(row_index)
                columns.append(column)
                values.append(float(value))
    matrix = coo_matrix(
        (values, (row_indices, columns)),
        shape=(len(equalities), variable_count),
    ).tocsr()
    bounds = [
        (0.0, 0.0) if variable in forbidden else (0.0, None)
        for variable in range(variable_count)
    ]
    result = linprog(
        np.zeros(variable_count),
        A_eq=matrix,
        b_eq=np.array([float(value) for value in rhs]),
        bounds=bounds,
        method="highs-ds",
    )
    if not result.success:
        return False, None, "float_infeasible"
    point = [F(value).limit_denominator(1_000_000) for value in result.x]
    exact = all(value >= 0 for value in point) and all(
        sum((coefficient * point[column] for column, coefficient in enumerate(row)), F(0))
        == target
        for row, target in zip(equalities, rhs, strict=True)
    ) and all(point[variable] == 0 for variable in forbidden)
    if exact:
        return True, point, "float_solution_exactly_reconstructed"

    active = [
        variable
        for variable, value in enumerate(result.x)
        if value > 1e-9 and variable not in forbidden
    ]
    reduced_equalities = [
        [row[variable] for variable in active] for row in equalities
    ]
    nonnegative_rows = []
    for variable in range(len(active)):
        row = [F(0)] * len(active)
        row[variable] = F(-1)
        nonnegative_rows.append(row)
    status, _, reduced_point = exact_linprog(
        [F(0)] * len(active),
        nonnegative_rows,
        [F(0)] * len(active),
        reduced_equalities,
        rhs,
    )
    if status != "optimal" or reduced_point is None:
        return False, None, "active_column_exact_repair_failed"
    point = [F(0)] * variable_count
    for variable, value in zip(active, reduced_point, strict=True):
        point[variable] = value
    return True, point, "active_column_exact_repair"


def terminals(payload: dict[str, Any]) -> tuple[int, list[dict[str, Any]]]:
    n = int(payload["n"])
    game_count = len(payload["family"]["games"])
    rows = [
        active
        for active in payload["exact_margin_dual_certificate"]["active_inequalities"]
        if active["row"][0] == "monotonicity"
    ]
    scale = lcm(*(F(active["weight"]).denominator for active in rows))
    divergence = [[F(0)] * n for _ in range(game_count)]
    for active in rows:
        _, lower, upper, _, player = active["row"]
        flow = -F(active["weight"])
        divergence[int(lower)][int(player)] += flow
        divergence[int(upper)][int(player)] -= flow
    result = []
    for node, vector in enumerate(divergence):
        if not any(vector):
            continue
        nonzero = {value for value in vector if value}
        if len(nonzero) != 1:
            raise RuntimeError(f"node {node} has mixed terminal divergence")
        coefficient = next(iter(nonzero))
        coalition = sum(
            1 << player for player, value in enumerate(vector) if value
        )
        result.append(
            {
                "node": node,
                "coefficient": coefficient * scale,
                "coalition": coalition,
            }
        )
    return n, result


def protected_reachability(
    payload: dict[str, Any], n: int, sources: list[dict[str, Any]]
) -> dict[tuple[int, int], set[int]]:
    adjacency = [[[] for _ in payload["family"]["games"]] for _ in range(n)]
    for edge in payload["family"]["edges"]:
        coalition = int(edge["coalition"])
        for player in range(n):
            if coalition >> player & 1:
                adjacency[player][int(edge["lower"])].append(int(edge["upper"]))
    reachable: dict[tuple[int, int], set[int]] = {}
    for source in sources:
        source_node = int(source["node"])
        for player in range(n):
            seen = {source_node}
            stack = [source_node]
            while stack:
                node = stack.pop()
                for successor in adjacency[player][node]:
                    if successor not in seen:
                        seen.add(successor)
                        stack.append(successor)
            reachable[source_node, player] = seen
    return reachable


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.archive.read_text())
    n, raw_terminals = terminals(payload)
    sources = [terminal for terminal in raw_terminals if terminal["coefficient"] > 0]
    sinks = [terminal for terminal in raw_terminals if terminal["coefficient"] < 0]
    for sink in sinks:
        sink["coefficient"] = -sink["coefficient"]

    reachable = protected_reachability(payload, n, sources)

    variable_count = len(sources) * len(sinks)
    equalities: list[list[F]] = []
    rhs: list[F] = []
    for source_index, source in enumerate(sources):
        for player in range(n):
            row = [F(0)] * variable_count
            for sink_index, sink in enumerate(sinks):
                if sink["coalition"] >> player & 1:
                    row[source_index * len(sinks) + sink_index] = F(1)
            equalities.append(row)
            rhs.append(
                F(source["coefficient"])
                if source["coalition"] >> player & 1
                else F(0)
            )
    for sink_index, sink in enumerate(sinks):
        row = [F(0)] * variable_count
        for source_index in range(len(sources)):
            row[source_index * len(sinks) + sink_index] = F(1)
        equalities.append(row)
        rhs.append(F(sink["coefficient"]))
    fork_algebraically_feasible = grouped_algebraic_solve(
        sources, sinks, n, "fork"
    )
    fork_algebraic_method = "exact_grouped_coalition_types"
    fork_incompatible_pairs = []
    fork_forbidden: set[int] = set()
    if fork_algebraically_feasible:
        for source_index, source in enumerate(sources):
            for sink_index, sink in enumerate(sinks):
                compatible = all(
                    int(sink["node"]) in reachable[int(source["node"]), player]
                    for player in range(n)
                    if int(sink["coalition"]) >> player & 1
                )
                if not compatible:
                    fork_forbidden.add(source_index * len(sinks) + sink_index)
                    fork_incompatible_pairs.append(
                        {"source": source["node"], "sink": sink["node"]}
                    )
        fork_feasible, point, fork_method = transport_solve(
            variable_count,
            equalities,
            rhs,
            fork_forbidden,
        )
    else:
        fork_feasible, point = False, None
        fork_method = "skipped_exact_grouped_infeasible"
    fork_decomposition = []
    if fork_feasible:
        assert point is not None
        for source_index, source in enumerate(sources):
            for sink_index, sink in enumerate(sinks):
                weight = point[source_index * len(sinks) + sink_index]
                if weight:
                    fork_decomposition.append(
                        {
                            "source": source["node"],
                            "sink": sink["node"],
                            "weight": str(weight),
                        }
                    )

    join_equalities: list[list[F]] = []
    join_rhs: list[F] = []
    for source_index, source in enumerate(sources):
        row = [F(0)] * variable_count
        for sink_index in range(len(sinks)):
            row[source_index * len(sinks) + sink_index] = F(1)
        join_equalities.append(row)
        join_rhs.append(F(source["coefficient"]))
    for sink_index, sink in enumerate(sinks):
        for player in range(n):
            row = [F(0)] * variable_count
            for source_index, source in enumerate(sources):
                if int(source["coalition"]) >> player & 1:
                    row[source_index * len(sinks) + sink_index] = F(1)
            join_equalities.append(row)
            join_rhs.append(
                F(sink["coefficient"])
                if int(sink["coalition"]) >> player & 1
                else F(0)
            )
    join_algebraically_feasible = grouped_algebraic_solve(
        sources, sinks, n, "join"
    )
    join_algebraic_method = "exact_grouped_coalition_types"
    join_incompatible_pairs = []
    join_forbidden: set[int] = set()
    if join_algebraically_feasible:
        for source_index, source in enumerate(sources):
            for sink_index, sink in enumerate(sinks):
                compatible = all(
                    int(sink["node"]) in reachable[int(source["node"]), player]
                    for player in range(n)
                    if int(source["coalition"]) >> player & 1
                )
                if not compatible:
                    join_forbidden.add(source_index * len(sinks) + sink_index)
                    join_incompatible_pairs.append(
                        {"source": source["node"], "sink": sink["node"]}
                    )
        join_feasible, join_point, join_method = transport_solve(
            variable_count,
            join_equalities,
            join_rhs,
            join_forbidden,
        )
    else:
        join_feasible, join_point = False, None
        join_method = "skipped_exact_grouped_infeasible"
    join_decomposition = []
    if join_feasible:
        assert join_point is not None
        for source_index, source in enumerate(sources):
            for sink_index, sink in enumerate(sinks):
                weight = join_point[source_index * len(sinks) + sink_index]
                if weight:
                    join_decomposition.append(
                        {
                            "source": source["node"],
                            "sink": sink["node"],
                            "weight": str(weight),
                        }
                    )
    feasible = fork_feasible or join_feasible
    result = {
        "status": "indicator_block_transport_audit",
        "source": str(args.archive),
        "terminal_count": len(raw_terminals),
        "source_count": len(sources),
        "sink_count": len(sinks),
        "algebraically_decomposable": (
            fork_algebraically_feasible or join_algebraically_feasible
        ),
        "block_transport_decomposable": feasible,
        "theorem_applies": feasible,
        "reachability_restricted": True,
        "fork_transport": {
            "algebraically_decomposable": fork_algebraically_feasible,
            "decomposable": fork_feasible,
            "algebraic_solve_method": fork_algebraic_method,
            "restricted_solve_method": fork_method,
            "incompatible_pairs": fork_incompatible_pairs,
            "decomposition": fork_decomposition,
        },
        "join_transport": {
            "algebraically_decomposable": join_algebraically_feasible,
            "decomposable": join_feasible,
            "algebraic_solve_method": join_algebraic_method,
            "restricted_solve_method": join_method,
            "incompatible_pairs": join_incompatible_pairs,
            "decomposition": join_decomposition,
        },
        "sources": [
            {
                "node": source["node"],
                "coalition": source["coalition"],
                "coefficient": str(source["coefficient"]),
            }
            for source in sources
        ],
        "sinks": [
            {
                "node": sink["node"],
                "coalition": sink["coalition"],
                "coefficient": str(sink["coefficient"]),
            }
            for sink in sinks
        ],
        "decomposition": (
            {"orientation": "fork", "weights": fork_decomposition}
            if fork_feasible
            else (
                {"orientation": "join", "weights": join_decomposition}
                if join_feasible
                else None
            )
        ),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
