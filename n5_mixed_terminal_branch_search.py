#!/usr/bin/env python3
"""Search all extreme lower-expectation branches at three-level terminals."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, csr_matrix, hstack, vstack

from n5_facet_search import load_facets


F = Fraction


def sparse_matrix(
    rows: list[dict[int, float]], column_count: int
) -> csr_matrix:
    row_indices = []
    columns = []
    values = []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                row_indices.append(row_index)
                columns.append(column)
                values.append(value)
    return coo_matrix(
        (values, (row_indices, columns)),
        shape=(len(rows), column_count),
    ).tocsr()


def extreme_representations(
    divergence: list[F], n: int
) -> list[tuple[F, tuple[tuple[int, F], ...]]]:
    proper = range(1, (1 << n) - 1)
    representations = set()
    for coalitions in itertools.combinations(proper, n - 1):
        matrix = np.asarray(
            [
                [float((coalition >> player) & 1) for coalition in coalitions]
                + [1.0]
                for player in range(n)
            ]
        )
        if abs(float(np.linalg.det(matrix))) < 0.5:
            continue
        solution = np.linalg.solve(
            matrix, np.asarray([float(value) for value in divergence])
        )
        rational = tuple(
            F(float(value)).limit_denominator(10_000) for value in solution
        )
        if any(value < 0 for value in rational[:-1]):
            continue
        if any(
            sum(
                (
                    rational[index]
                    * ((coalition >> player) & 1)
                    for index, coalition in enumerate(coalitions)
                ),
                rational[-1],
            )
            != divergence[player]
            for player in range(n)
        ):
            continue
        representations.add(
            (
                rational[-1],
                tuple(
                    (coalition, rational[index])
                    for index, coalition in enumerate(coalitions)
                    if rational[index]
                ),
            )
        )
    return sorted(representations, key=str)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument(
        "--omit-mixed-node", type=int, action="append", default=[]
    )
    parser.add_argument(
        "--branch-encoding",
        choices=("one-hot", "binary", "hybrid", "perspective"),
        default="one-hot",
    )
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    node_index = {node: index for index, node in enumerate(node_ids)}
    worth_variable_count = len(node_ids) * coalition_count

    def column(node: int, coalition: int) -> int:
        return node_index[node] * coalition_count + coalition

    inequalities: list[dict[int, float]] = []
    inequality_rhs: list[float] = []
    equalities: list[dict[int, float]] = []
    equality_rhs: list[float] = []
    facets, _ = load_facets()
    for node in node_ids:
        equalities.append({column(node, 0): 1.0})
        equality_rhs.append(0.0)
        equalities.append({column(node, grand): 1.0})
        equality_rhs.append(1.0)
        for facet in facets:
            inequalities.append(
                {
                    column(node, coalition): -float(coefficient)
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                }
            )
            inequality_rhs.append(0.0)
        for coalition in range(coalition_count):
            for player in range(n):
                if coalition >> player & 1:
                    continue
                successor = coalition | (1 << player)
                inequalities.append(
                    {
                        column(node, coalition): 1.0,
                        column(node, successor): -1.0,
                    }
                )
                inequality_rhs.append(0.0)

    protected_players: dict[tuple[int, int], set[int]] = defaultdict(set)
    for path in report["terminal_path_decomposition"]:
        protected_players[
            int(path["source"]), int(path["sink"])
        ].add(int(path["player"]))
    for (source, sink), players in protected_players.items():
        player_mask = sum(1 << player for player in players)
        for coalition in range(coalition_count):
            row = {
                column(source, coalition): 1.0,
                column(sink, coalition): -1.0,
            }
            if coalition & player_mask != player_mask:
                equalities.append(row)
                equality_rhs.append(0.0)
            else:
                inequalities.append(row)
                inequality_rhs.append(0.0)

    all_mixed = [
        terminal for terminal in terminals if terminal.get("coefficient") is None
    ]
    omitted_mixed_nodes = set(args.omit_mixed_node)
    mixed = [
        terminal
        for terminal in all_mixed
        if int(terminal["node"]) not in omitted_mixed_nodes
    ]
    if not mixed:
        raise RuntimeError("terminal report has no non-affine terminals")
    representations = {
        int(terminal["node"]): extreme_representations(
            [F(value) for value in terminal["scaled_divergence"]], n
        )
        for terminal in mixed
    }
    if any(not choices for choices in representations.values()):
        raise RuntimeError("a mixed terminal has no extreme representation")

    z_offset = worth_variable_count
    auxiliary_offset = z_offset + len(mixed)
    choice_ranges = {}
    bit_ranges = {}
    copy_offsets = {}
    next_variable = auxiliary_offset
    for terminal in mixed:
        node = int(terminal["node"])
        bit_count = math.ceil(math.log2(len(representations[node])))
        if args.branch_encoding in ("one-hot", "hybrid", "perspective"):
            choice_ranges[node] = (
                next_variable,
                next_variable + len(representations[node]),
            )
            next_variable += len(representations[node])
        if args.branch_encoding in ("binary", "hybrid", "perspective"):
            bit_ranges[node] = (
                next_variable,
                next_variable + bit_count,
            )
            next_variable += bit_count
        if args.branch_encoding == "perspective":
            copy_offsets[node] = next_variable
            next_variable += len(representations[node]) * coalition_count
    total_variables = next_variable
    branch_rows: list[dict[int, float]] = []
    branch_rhs: list[float] = []
    selection_rows: list[dict[int, float]] = []
    selection_lower: list[float] = []
    selection_upper: list[float] = []
    for mixed_index, terminal in enumerate(mixed):
        node = int(terminal["node"])
        choice_start, choice_stop = choice_ranges.get(node, (0, 0))
        bit_start, bit_stop = bit_ranges.get(node, (0, 0))
        if args.branch_encoding == "perspective":
            choices = representations[node]
            copy_offset = copy_offsets[node]

            selection_rows.append(
                {
                    choice: 1.0
                    for choice in range(choice_start, choice_stop)
                }
            )
            selection_lower.append(1.0)
            selection_upper.append(1.0)
            for bit, binary in enumerate(range(bit_start, bit_stop)):
                row = {binary: 1.0}
                for choice_index, choice in enumerate(
                    range(choice_start, choice_stop)
                ):
                    if choice_index >> bit & 1:
                        row[choice] = -1.0
                selection_rows.append(row)
                selection_lower.append(0.0)
                selection_upper.append(0.0)

            objective_row = {z_offset + mixed_index: 1.0}
            for choice_index, (mu, coalitions) in enumerate(choices):
                choice = choice_start + choice_index
                objective_row[choice] = -float(mu)
                local_offset = copy_offset + choice_index * coalition_count
                for coalition, weight in coalitions:
                    objective_row[local_offset + coalition] = -float(weight)
            selection_rows.append(objective_row)
            selection_lower.append(0.0)
            selection_upper.append(0.0)

            for coalition in range(coalition_count):
                row = {column(node, coalition): 1.0}
                for choice_index in range(len(choices)):
                    local_offset = (
                        copy_offset + choice_index * coalition_count
                    )
                    row[local_offset + coalition] = -1.0
                selection_rows.append(row)
                selection_lower.append(0.0)
                selection_upper.append(0.0)

            for choice_index in range(len(choices)):
                choice = choice_start + choice_index
                local_offset = copy_offset + choice_index * coalition_count
                selection_rows.append({local_offset: 1.0})
                selection_lower.append(0.0)
                selection_upper.append(0.0)
                selection_rows.append(
                    {
                        local_offset + grand: 1.0,
                        choice: -1.0,
                    }
                )
                selection_lower.append(0.0)
                selection_upper.append(0.0)
                for facet in facets:
                    branch_rows.append(
                        {
                            local_offset + coalition: -float(coefficient)
                            for coalition, coefficient in enumerate(facet)
                            if coefficient
                        }
                    )
                    branch_rhs.append(0.0)
                for coalition in range(coalition_count):
                    for player in range(n):
                        if coalition >> player & 1:
                            continue
                        successor = coalition | (1 << player)
                        branch_rows.append(
                            {
                                local_offset + coalition: 1.0,
                                local_offset + successor: -1.0,
                            }
                        )
                        branch_rhs.append(0.0)
            continue

        one_hot = {}
        for choice_index, (mu, coalitions) in enumerate(
            representations[node]
        ):
            big_m = -float(mu)
            if big_m <= 0:
                raise RuntimeError("mixed branch has no valid grand offset")
            row = {
                z_offset + mixed_index: 1.0,
                column(node, grand): -float(mu),
            }
            for coalition, weight in coalitions:
                row[column(node, coalition)] = (
                    row.get(column(node, coalition), 0.0) - float(weight)
                )
            if args.branch_encoding in ("one-hot", "hybrid"):
                choice = choice_start + choice_index
                one_hot[choice] = 1.0
                row[choice] = big_m
                rhs = big_m
            else:
                one_bits = 0
                for bit, binary in enumerate(range(bit_start, bit_stop)):
                    if choice_index >> bit & 1:
                        row[binary] = big_m
                        one_bits += 1
                    else:
                        row[binary] = -big_m
                rhs = big_m * one_bits
            branch_rows.append(row)
            branch_rhs.append(rhs)
        if args.branch_encoding in ("one-hot", "hybrid"):
            selection_rows.append(one_hot)
            selection_lower.append(1.0)
            selection_upper.append(1.0)
        if args.branch_encoding == "binary":
            selection_rows.append(
                {
                    binary: float(1 << bit)
                    for bit, binary in enumerate(range(bit_start, bit_stop))
                }
            )
            selection_lower.append(-np.inf)
            selection_upper.append(float(len(representations[node]) - 1))
        elif args.branch_encoding == "hybrid":
            for bit, binary in enumerate(range(bit_start, bit_stop)):
                row = {binary: 1.0}
                for choice_index, choice in enumerate(
                    range(choice_start, choice_stop)
                ):
                    if choice_index >> bit & 1:
                        row[choice] = -1.0
                selection_rows.append(row)
                selection_lower.append(0.0)
                selection_upper.append(0.0)

    pad_columns = total_variables - worth_variable_count
    base_ub = hstack(
        (
            sparse_matrix(inequalities, worth_variable_count),
            csr_matrix((len(inequalities), pad_columns)),
        ),
        format="csr",
    )
    base_eq = hstack(
        (
            sparse_matrix(equalities, worth_variable_count),
            csr_matrix((len(equalities), pad_columns)),
        ),
        format="csr",
    )
    matrix = vstack(
        (
            base_ub,
            sparse_matrix(branch_rows, total_variables),
            base_eq,
            sparse_matrix(selection_rows, total_variables),
        ),
        format="csr",
    )
    lower = np.concatenate(
        (
            np.full(len(inequalities), -np.inf),
            np.full(len(branch_rows), -np.inf),
            np.asarray(equality_rhs),
            np.asarray(selection_lower),
        )
    )
    upper = np.concatenate(
        (
            np.asarray(inequality_rhs),
            np.asarray(branch_rhs),
            np.asarray(equality_rhs),
            np.asarray(selection_upper),
        )
    )

    objective = np.zeros(total_variables)
    objective_constant = 0.0
    mixed_nodes = {int(terminal["node"]) for terminal in all_mixed}
    for terminal in terminals:
        node = int(terminal["node"])
        if node in mixed_nodes:
            continue
        coefficient = F(terminal["coefficient"])
        coalition = int(terminal["coalition"])
        if coefficient > 0:
            objective[column(node, coalition)] -= float(coefficient)
        else:
            mass = -coefficient
            objective[column(node, grand ^ coalition)] -= float(mass)
            objective_constant -= float(mass)
    for mixed_index in range(len(mixed)):
        objective[z_offset + mixed_index] = -1.0

    lower_bounds = np.full(total_variables, -np.inf)
    upper_bounds = np.full(total_variables, np.inf)
    lower_bounds[z_offset:auxiliary_offset] = -2.0
    upper_bounds[z_offset:auxiliary_offset] = 0.0
    lower_bounds[auxiliary_offset:] = 0.0
    upper_bounds[auxiliary_offset:] = 1.0
    integrality = np.zeros(total_variables)
    integer_variable_count = 0
    integer_ranges = (
        choice_ranges.values()
        if args.branch_encoding == "one-hot"
        else bit_ranges.values()
    )
    for start, stop in integer_ranges:
        integrality[start:stop] = 1
        integer_variable_count += stop - start
    result = milp(
        objective,
        integrality=integrality,
        bounds=Bounds(lower_bounds, upper_bounds),
        constraints=LinearConstraint(matrix, lower, upper),
        options={"time_limit": args.time_limit, "mip_rel_gap": 0.0},
    )
    if result.x is None:
        raise RuntimeError(result.message)

    selected = {}
    for terminal in mixed:
        node = int(terminal["node"])
        if args.branch_encoding in ("one-hot", "hybrid", "perspective"):
            start, stop = choice_ranges[node]
            choice = int(np.argmax(result.x[start:stop]))
        else:
            start, stop = bit_ranges[node]
            choice = sum(
                (1 << bit) * int(round(result.x[binary]))
                for bit, binary in enumerate(range(start, stop))
            )
        mu, coalitions = representations[node][choice]
        selected[str(node)] = {
            "choice": choice,
            "mu": str(mu),
            "coalitions": [
                {"coalition": coalition, "weight": str(weight)}
                for coalition, weight in coalitions
            ],
        }
    payload: dict[str, Any] = {
        "status": "mixed_terminal_branch_milp",
        "source": str(args.terminal_report),
        "solver_status": int(result.status),
        "solver_message": result.message,
        "success": bool(result.success),
        "maximum_lower_expectation_objective_float": (
            -float(result.fun) + objective_constant
        ),
        "mip_gap": (
            float(result.mip_gap) if result.mip_gap is not None else None
        ),
        "mip_node_count": (
            int(result.mip_node_count)
            if result.mip_node_count is not None
            else None
        ),
        "mixed_terminal_count": len(mixed),
        "omitted_mixed_nodes": sorted(omitted_mixed_nodes),
        "representation_counts": {
            str(node): len(choices)
            for node, choices in representations.items()
        },
        "selected_representations": selected,
        "big_m_rule": "branch-specific -mu",
        "branch_encoding": args.branch_encoding,
        "binary_variable_count": integer_variable_count,
        "continuous_choice_variable_count": (
            sum(stop - start for start, stop in choice_ranges.values())
            if args.branch_encoding in ("hybrid", "perspective")
            else 0
        ),
        "perspective_copy_variable_count": (
            sum(
                len(representations[node]) * coalition_count
                for node in copy_offsets
            )
        ),
        "interpretation": (
            "Each mixed terminal lower expectation is the maximum over its "
            "extreme nonnegative core-dual representations. The MILP selects "
            "one branch per terminal while optimizing all terminal games."
        ),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
