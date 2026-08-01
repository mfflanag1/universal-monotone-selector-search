#!/usr/bin/env python3
"""Exactify a perspective certificate after eliminating local block weights."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import qr
from scipy.sparse import coo_matrix
from sympy import Matrix, Rational, SparseMatrix, linsolve

from n5_facet_search import load_facets
from n5_mixed_terminal_branch_search import extreme_representations


F = Fraction


def exact_square_solve(
    equations: list[dict[int, F]],
    rhs: list[F],
    variable_count: int,
    reference: list[F] | None = None,
) -> tuple[list[F], dict[str, Any]]:
    row_indices = []
    columns = []
    values = []
    for row_index, row in enumerate(equations):
        for column, value in row.items():
            if value:
                row_indices.append(row_index)
                columns.append(column)
                values.append(float(value))
    matrix = coo_matrix(
        (values, (row_indices, columns)),
        shape=(len(equations), variable_count),
    ).tocsr()
    _, triangular, pivots = qr(
        matrix.transpose().toarray(),
        mode="economic",
        pivoting=True,
        check_finite=False,
    )
    diagonal = np.abs(np.diag(triangular))
    tolerance = (
        max(matrix.shape)
        * np.finfo(float).eps
        * (diagonal[0] if len(diagonal) else 0.0)
    )
    rank = int(np.count_nonzero(diagonal > tolerance))
    diagnostics = {
        "equation_count": len(equations),
        "variable_count": variable_count,
        "numeric_rank": rank,
        "rank_tolerance": tolerance,
        "smallest_selected_diagonal": (
            float(diagonal[rank - 1]) if rank else None
        ),
    }
    if rank < variable_count and reference is not None:
        _, column_triangular, column_pivots = qr(
            matrix.toarray(),
            mode="economic",
            pivoting=True,
            check_finite=False,
        )
        column_diagonal = np.abs(np.diag(column_triangular))
        column_rank = int(
            np.count_nonzero(column_diagonal > tolerance)
        )
        if column_rank != rank:
            raise RuntimeError("row and column rank diagnostics disagree")
        fixed_variables = [
            int(column) for column in column_pivots[rank:]
        ]
        augmented_equations = equations + [
            {column: F(1)} for column in fixed_variables
        ]
        augmented_rhs = rhs + [
            reference[column] for column in fixed_variables
        ]
        solution, augmented_diagnostics = exact_square_solve(
            augmented_equations,
            augmented_rhs,
            variable_count,
            None,
        )
        augmented_diagnostics["fixed_free_variables"] = fixed_variables
        augmented_diagnostics["fixed_free_values"] = [
            str(reference[column]) for column in fixed_variables
        ]
        return solution, augmented_diagnostics
    if rank != variable_count:
        raise RuntimeError(
            f"reduced system rank {rank} != variable count {variable_count}"
        )
    selected = [int(value) for value in pivots[:rank]]
    square = matrix[selected, :].tocoo()
    exact_matrix = SparseMatrix(
        rank,
        rank,
        {
            (int(row), int(column)): Rational(
                equations[selected[int(row)]][int(column)].numerator,
                equations[selected[int(row)]][int(column)].denominator,
            )
            for row, column in zip(square.row, square.col, strict=True)
        },
    )
    exact_rhs = Matrix(
        [
            Rational(rhs[row].numerator, rhs[row].denominator)
            for row in selected
        ]
    )
    solution = next(iter(linsolve((exact_matrix, exact_rhs))))
    weights = [F(int(value.p), int(value.q)) for value in solution]
    for row, expected in zip(equations, rhs, strict=True):
        actual = sum(
            (coefficient * weights[column] for column, coefficient in row.items()),
            F(0),
        )
        if actual != expected:
            raise RuntimeError("reduced exact solution fails a full equation")
    return weights, diagnostics


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("float_support", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--diagnostics-only", action="store_true")
    parser.add_argument("--anchor-search", type=Path)
    parser.add_argument("--simplex-reference-value", type=F)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    support = json.loads(args.float_support.read_text())
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    node_index = {node: index for index, node in enumerate(node_ids)}
    base_variable_count = len(node_ids) * coalition_count
    mixed_terminals = [
        terminal
        for terminal in terminals
        if terminal.get("coefficient") is None
        and len(set(terminal["scaled_divergence"])) > 2
    ]
    mixed_nodes = sorted(int(terminal["node"]) for terminal in mixed_terminals)
    facets, _ = load_facets()
    representations = {
        int(terminal["node"]): extreme_representations(
            [F(value) for value in terminal["scaled_divergence"]], n
        )
        for terminal in mixed_terminals
    }

    def base_column(node: int, coalition: int) -> int:
        return node_index[node] * coalition_count + coalition

    def global_row(metadata: list[Any]) -> tuple[dict[int, F], F]:
        row_type = metadata[0]
        if row_type == "exact_cone":
            _, node, facet_index = metadata
            facet = facets[int(facet_index)]
            return (
                {
                    base_column(int(node), coalition): -F(coefficient)
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                },
                F(0),
            )
        if row_type == "game_monotonicity":
            _, node, coalition, successor = metadata
            return (
                {
                    base_column(int(node), int(coalition)): F(1),
                    base_column(int(node), int(successor)): F(-1),
                },
                F(0),
            )
        if row_type == "directed_monotonicity":
            _, source, sink, coalition = metadata
            return (
                {
                    base_column(int(source), int(coalition)): F(1),
                    base_column(int(sink), int(coalition)): F(-1),
                },
                F(0),
            )
        if row_type == "empty":
            _, node = metadata
            return ({base_column(int(node), 0): F(1)}, F(0))
        if row_type == "grand":
            _, node = metadata
            return ({base_column(int(node), grand): F(1)}, F(1))
        if row_type == "protected_invariance":
            _, source, sink, coalition = metadata
            return (
                {
                    base_column(int(source), int(coalition)): F(1),
                    base_column(int(sink), int(coalition)): F(-1),
                },
                F(0),
            )
        raise RuntimeError(f"unsupported global row {row_type}")

    global_inequalities = [
        active
        for active in support["inequalities"]
        if not active["row"][0].startswith("perspective_")
    ]
    global_equalities = [
        active
        for active in support["equalities"]
        if not active["row"][0].startswith("perspective_")
    ]
    global_entries = global_inequalities + global_equalities
    global_inequality_count = len(global_inequalities)
    global_rows = [global_row(active["row"])[0] for active in global_entries]
    global_rhs = [global_row(active["row"])[1] for active in global_entries]

    local_support_rows: dict[tuple[int, int], list[dict[str, Any]]] = (
        defaultdict(list)
    )
    for active in support["inequalities"]:
        metadata = active["row"]
        if metadata[0] == "perspective_copy_exact_cone":
            _, node, choice, facet_index = metadata
            local_support_rows[int(node), int(choice)].append(
                {
                    "row": metadata,
                    "inequality": True,
                    "vector": [
                        F(-coefficient)
                        for coefficient in facets[int(facet_index)]
                    ]
                    + [F(0)],
                }
            )
        elif metadata[0] == "perspective_copy_game_monotonicity":
            _, node, choice, coalition, successor = metadata
            vector = [F(0)] * (coalition_count + 1)
            vector[int(coalition)] = F(1)
            vector[int(successor)] = F(-1)
            local_support_rows[int(node), int(choice)].append(
                {
                    "row": metadata,
                    "inequality": True,
                    "vector": vector,
                }
            )
        elif metadata[0] == "perspective_choice_nonnegative":
            _, node, choice = metadata
            vector = [F(0)] * (coalition_count + 1)
            vector[-1] = F(-1)
            local_support_rows[int(node), int(choice)].append(
                {
                    "row": metadata,
                    "inequality": True,
                    "vector": vector,
                }
            )
        elif metadata[0].startswith("perspective_"):
            raise RuntimeError(f"unsupported local inequality {metadata[0]}")
    for active in support["equalities"]:
        metadata = active["row"]
        if metadata[0] == "perspective_copy_empty":
            _, node, choice = metadata
            vector = [F(0)] * (coalition_count + 1)
            vector[0] = F(1)
            local_support_rows[int(node), int(choice)].append(
                {
                    "row": metadata,
                    "inequality": False,
                    "vector": vector,
                }
            )
        elif metadata[0] == "perspective_copy_grand":
            _, node, choice = metadata
            vector = [F(0)] * (coalition_count + 1)
            vector[grand] = F(1)
            vector[-1] = F(-1)
            local_support_rows[int(node), int(choice)].append(
                {
                    "row": metadata,
                    "inequality": False,
                    "vector": vector,
                }
            )
        elif metadata[0] in (
            "perspective_link",
            "perspective_simplex",
        ):
            continue
        elif metadata[0].startswith("perspective_"):
            raise RuntimeError(f"unsupported local equality {metadata[0]}")

    global_offset = 0
    q_offsets = {
        node: len(global_entries) + index * coalition_count
        for index, node in enumerate(mixed_nodes)
    }
    sigma_offset = len(global_entries) + len(mixed_nodes) * coalition_count
    sigma_columns = {
        node: sigma_offset + index for index, node in enumerate(mixed_nodes)
    }
    reduced_variable_count = sigma_offset + len(mixed_nodes)

    objective = [F(0)] * base_variable_count
    for terminal in terminals:
        node = int(terminal["node"])
        if node in mixed_nodes:
            continue
        raw_coefficient = terminal.get("coefficient")
        if raw_coefficient is not None:
            coefficient = F(raw_coefficient)
            coalition = int(terminal["coalition"])
            if coefficient > 0:
                objective[base_column(node, coalition)] -= coefficient
            elif coefficient < 0:
                objective[
                    base_column(node, grand ^ coalition)
                ] -= -coefficient
        else:
            divergence = [F(value) for value in terminal["scaled_divergence"]]
            low, high = sorted(set(divergence))
            coalition = sum(
                1 << player
                for player, value in enumerate(divergence)
                if value == high
            )
            objective[base_column(node, coalition)] -= high - low
            objective[base_column(node, grand)] -= low

    equations: list[dict[int, F]] = [
        {} for _ in range(base_variable_count)
    ]
    rhs = list(objective)
    for support_column, row in enumerate(global_rows):
        for variable, coefficient in row.items():
            equations[variable][global_offset + support_column] = coefficient
    for node in mixed_nodes:
        for coalition in range(coalition_count):
            equations[base_column(node, coalition)][
                q_offsets[node] + coalition
            ] = F(1)
    objective_equation: dict[int, F] = {}
    for support_column, bound in enumerate(global_rhs):
        if bound:
            objective_equation[global_offset + support_column] = bound
    for node in mixed_nodes:
        objective_equation[sigma_columns[node]] = F(1)
    equations.append(objective_equation)
    objective_constant = sum(
        (
            F(terminal["coefficient"])
            for terminal in terminals
            if terminal.get("coefficient") is not None
            and F(terminal["coefficient"]) < 0
        ),
        F(0),
    )
    rhs.append(objective_constant)

    for node in mixed_nodes:
        for choice, (mu, coalitions) in enumerate(representations[node]):
            metadata = local_support_rows[node, choice]
            local_matrix_columns = [
                entry["vector"] for entry in metadata
            ]
            local_matrix = Matrix(
                [
                    [
                        Rational(
                            local_matrix_columns[column][row].numerator,
                            local_matrix_columns[column][row].denominator,
                        )
                        for column in range(len(local_matrix_columns))
                    ]
                    for row in range(coalition_count + 1)
                ]
            )
            left_nullspace = local_matrix.transpose().nullspace()
            weights = [F(0)] * coalition_count
            for coalition, weight in coalitions:
                weights[coalition] = weight
            base_target = [-weight for weight in weights] + [-mu]
            for null_vector in left_nullspace:
                equation: dict[int, F] = {}
                for coalition in range(coalition_count):
                    coefficient = F(
                        int(null_vector[coalition].p),
                        int(null_vector[coalition].q),
                    )
                    if coefficient:
                        equation[q_offsets[node] + coalition] = coefficient
                sigma_coefficient = -F(
                    int(null_vector[-1].p), int(null_vector[-1].q)
                )
                if sigma_coefficient:
                    equation[sigma_columns[node]] = sigma_coefficient
                expected = -sum(
                    (
                        F(int(value.p), int(value.q)) * target
                        for value, target in zip(
                            null_vector, base_target, strict=True
                        )
                    ),
                    F(0),
                )
                equations.append(equation)
                rhs.append(expected)

    if args.anchor_search is not None:
        anchor = json.loads(args.anchor_search.read_text())[
            "best_complete_anchor"
        ]
        fixed_node = min(mixed_nodes)
        equations.append({sigma_columns[fixed_node]: F(1)})
        rhs.append(F(anchor["simplex"][str(fixed_node)]))

    if args.diagnostics_only:
        print(
            json.dumps(
                {
                    "reduced_variable_count": reduced_variable_count,
                    "reduced_equation_count": len(equations),
                    "global_support_rows": len(global_entries),
                    "local_block_count": len(local_support_rows),
                },
                indent=2,
            )
        )
        return 0

    reduced_solution, diagnostics = exact_square_solve(
        equations,
        rhs,
        reduced_variable_count,
        [
            F(active["weight"]).limit_denominator(10_000_000)
            for active in global_entries
        ]
        + [
            next(
                (
                    F(active["weight"]).limit_denominator(10_000_000)
                    for active in support["equalities"]
                    if active["row"]
                    == ["perspective_link", node, coalition]
                ),
                F(0),
            )
            for node in mixed_nodes
            for coalition in range(coalition_count)
        ]
        + [
            (
                args.simplex_reference_value
                if args.simplex_reference_value is not None
                else next(
                    (
                        F(active["weight"]).limit_denominator(10_000_000)
                        for active in support["equalities"]
                        if active["row"] == ["perspective_simplex", node]
                    ),
                    F(0),
                )
            )
            for node in mixed_nodes
        ],
    )
    if any(
        reduced_solution[index] > 0
        for index in range(global_inequality_count)
    ):
        raise RuntimeError("reduced global inequality sign check failed")
    anchors = {
        node: [
            reduced_solution[q_offsets[node] + coalition]
            for coalition in range(coalition_count)
        ]
        for node in mixed_nodes
    }
    simplex = {
        node: reduced_solution[sigma_columns[node]]
        for node in mixed_nodes
    }

    local_certificates: dict[str, list[dict[str, Any]]] = {}
    local_inequality_count = 0
    local_equality_count = 0
    for node in mixed_nodes:
        node_certificates = []
        for choice, (mu, coalitions) in enumerate(representations[node]):
            metadata = local_support_rows[node, choice]
            local_columns = [entry["vector"] for entry in metadata]
            target = list(anchors[node])
            weights = [F(0)] * coalition_count
            for coalition, weight in coalitions:
                weights[coalition] = weight
            target = [
                anchor - weight
                for anchor, weight in zip(target, weights, strict=True)
            ]
            target.append(-mu - simplex[node])
            local_equations = [
                {
                    column: local_columns[column][row]
                    for column in range(len(local_columns))
                    if local_columns[column][row]
                }
                for row in range(coalition_count + 1)
            ]
            local_solution, _ = exact_square_solve(
                local_equations, target, len(local_columns)
            )
            for entry, value in zip(metadata, local_solution, strict=True):
                if entry["inequality"] and value > 0:
                    raise RuntimeError(
                        f"local inequality sign failed at {node}/{choice}"
                    )
            active_rows = [
                {
                    "row": entry["row"],
                    "weight": str(value),
                }
                for entry, value in zip(
                    metadata, local_solution, strict=True
                )
                if value
            ]
            local_inequality_count += sum(
                bool(entry["inequality"]) and bool(value)
                for entry, value in zip(
                    metadata, local_solution, strict=True
                )
            )
            local_equality_count += sum(
                not bool(entry["inequality"]) and bool(value)
                for entry, value in zip(
                    metadata, local_solution, strict=True
                )
            )
            node_certificates.append(
                {
                    "choice": choice,
                    "mu": str(mu),
                    "coalitions": [
                        {
                            "coalition": coalition,
                            "weight": str(weight),
                        }
                        for coalition, weight in coalitions
                    ],
                    "active_rows": active_rows,
                }
            )
        local_certificates[str(node)] = node_certificates

    global_weights = reduced_solution[: len(global_entries)]
    global_certificate = {
        "active_inequalities": [
            {
                "row": active["row"],
                "weight": str(weight),
            }
            for active, weight in zip(
                global_inequalities,
                global_weights[:global_inequality_count],
                strict=True,
            )
            if weight
        ],
        "active_equalities": [
            {
                "row": active["row"],
                "weight": str(weight),
            }
            for active, weight in zip(
                global_equalities,
                global_weights[global_inequality_count:],
                strict=True,
            )
            if weight
        ],
    }
    payload = {
        "status": "mixed_terminal_perspective_reduced_exact",
        "source": str(args.terminal_report),
        "float_support": str(args.float_support),
        "diagnostics": diagnostics,
        "objective_constant": str(objective_constant),
        "anchors": {
            str(node): [str(value) for value in values]
            for node, values in anchors.items()
        },
        "simplex": {
            str(node): str(value) for node, value in simplex.items()
        },
        "global_certificate": global_certificate,
        "local_block_count": len(local_support_rows),
        "local_inequality_row_count": local_inequality_count,
        "local_equality_row_count": local_equality_count,
        "local_certificates": local_certificates,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "reduced_variable_count": diagnostics["variable_count"],
                "reduced_equation_count": diagnostics["equation_count"],
                "global_inequality_rows": len(
                    global_certificate["active_inequalities"]
                ),
                "global_equality_rows": len(
                    global_certificate["active_equalities"]
                ),
                "local_block_count": payload["local_block_count"],
                "local_inequality_row_count": local_inequality_count,
                "local_equality_row_count": local_equality_count,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
