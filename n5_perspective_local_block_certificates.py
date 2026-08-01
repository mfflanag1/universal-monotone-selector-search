#!/usr/bin/env python3
"""Exactify every local representation block at rational perspective anchors."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import qr
from scipy.optimize import linprog
from sympy import Matrix, Rational, SparseMatrix, linsolve

from n5_facet_search import load_facets
from n5_mixed_terminal_branch_search import (
    extreme_representations,
    sparse_matrix,
)
from n5_terminal_topology_relaxation import reconstruct_exact_dual


F = Fraction


def exactify_basis(
    solve: object,
    matrix: object,
    target: list[F],
    inequality_count: int,
    inequality_metadata: list[tuple[object, ...]],
    equality_metadata: list[tuple[object, ...]],
) -> dict[str, object]:
    active = [
        index
        for index, value in enumerate(solve.x)
        if abs(float(value)) > 1e-12
    ]
    support_matrix = matrix[:, active].tocsr()
    _, triangular, pivots = qr(
        support_matrix.transpose().toarray(),
        mode="economic",
        pivoting=True,
        check_finite=False,
    )
    diagonal = np.abs(np.diag(triangular))
    tolerance = (
        max(support_matrix.shape)
        * np.finfo(float).eps
        * (diagonal[0] if len(diagonal) else 0.0)
    )
    rank = int(np.count_nonzero(diagonal > tolerance))
    if rank != len(active):
        raise RuntimeError("local active support is rank deficient")
    selected = [int(value) for value in pivots[:rank]]
    square = support_matrix[selected, :].tocoo()
    exact_square = SparseMatrix(
        rank,
        rank,
        {
            (int(row), int(column)): int(round(value))
            for row, column, value in zip(
                square.row, square.col, square.data, strict=True
            )
        },
    )
    exact_rhs = Matrix(
        [
            Rational(target[row].numerator, target[row].denominator)
            for row in selected
        ]
    )
    solution = next(iter(linsolve((exact_square, exact_rhs))))
    active_weights = [
        F(int(value.p), int(value.q)) for value in solution
    ]
    weights = [F(0)] * matrix.shape[1]
    for index, weight in zip(active, active_weights, strict=True):
        weights[index] = weight
    if any(weight > 0 for weight in weights[:inequality_count]):
        raise RuntimeError("local exact basis violates an inequality sign")
    exact_lhs = [F(0)] * matrix.shape[0]
    coo = matrix.tocoo()
    for row, column, coefficient in zip(
        coo.row, coo.col, coo.data, strict=True
    ):
        exact_lhs[int(row)] += weights[int(column)] * F(
            int(round(coefficient))
        )
    if exact_lhs != target:
        raise RuntimeError("local exact basis fails full stationarity")
    return {
        "reconstruction_method": "exact_active_basis_elimination",
        "variable_objective": "0",
        "lower_expectation_objective_with_constant": None,
        "active_inequalities": [
            {
                "row": list(inequality_metadata[index]),
                "weight": str(weight),
            }
            for index, weight in enumerate(weights[:inequality_count])
            if weight
        ],
        "active_equalities": [
            {
                "row": list(equality_metadata[index]),
                "weight": str(weight),
            }
            for index, weight in enumerate(weights[inequality_count:])
            if weight
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("global_support", type=Path)
    parser.add_argument("--simplex-283", default="2/3")
    parser.add_argument("--simplex-285", default="1")
    parser.add_argument("--anchor-search", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    global_payload = json.loads(args.global_support.read_text())
    links = {
        int(node): [F(value) for value in values]
        for node, values in global_payload["links"].items()
    }
    if args.anchor_search is None:
        simplex = {283: F(args.simplex_283), 285: F(args.simplex_285)}
    else:
        anchor = json.loads(args.anchor_search.read_text())[
            "best_complete_anchor"
        ]
        simplex = {
            int(node): F(value)
            for node, value in anchor["simplex"].items()
        }
        adjusted_node = max(simplex)
        simplex[adjusted_node] = (
            F(-80)
            - F(global_payload["global_variable_objective"])
            - sum(
                (
                    value
                    for node, value in simplex.items()
                    if node != adjusted_node
                ),
                F(0),
            )
        )
    if set(simplex) != set(links):
        raise RuntimeError("simplex nodes do not match global mixed nodes")
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    facets, _ = load_facets()
    representations = {
        int(terminal["node"]): extreme_representations(
            [F(value) for value in terminal["scaled_divergence"]], n
        )
        for terminal in terminals
        if int(terminal["node"]) in links
    }

    choice_variable = coalition_count
    variable_count = coalition_count + 1
    inequalities = [
        {
            coalition: -float(coefficient)
            for coalition, coefficient in enumerate(facet)
            if coefficient
        }
        for facet in facets
    ]
    inequality_metadata = [
        ("local_exact_cone", facet_index)
        for facet_index in range(len(facets))
    ]
    inequalities.append({choice_variable: -1.0})
    inequality_metadata.append(("local_choice_nonnegative",))
    inequality_rhs = [0.0] * len(inequalities)
    equalities = [
        {0: 1.0},
        {grand: 1.0, choice_variable: -1.0},
    ]
    equality_rhs = [0.0, 0.0]
    equality_metadata = [
        ("local_empty",),
        ("local_grand",),
    ]
    all_rows = inequalities + equalities
    matrix = sparse_matrix(all_rows, variable_count).transpose().tocsr()
    bounds = (
        [(None, 0.0)] * len(inequalities)
        + [(None, None)] * len(equalities)
    )

    certificates: dict[str, list[dict[str, object]]] = {}
    total_inequality_rows = 0
    total_equality_rows = 0
    for node, choices in representations.items():
        node_certificates = []
        for choice, (mu, coalitions) in enumerate(choices):
            target = [links[node][coalition] for coalition in range(coalition_count)]
            for coalition, weight in coalitions:
                target[coalition] -= weight
            target.append(-mu - simplex[node])
            solve = linprog(
                np.zeros(len(all_rows)),
                A_eq=matrix,
                b_eq=np.asarray([float(value) for value in target]),
                bounds=bounds,
                method="highs-ds",
                options={
                    "dual_feasibility_tolerance": 1e-10,
                    "primal_feasibility_tolerance": 1e-10,
                },
            )
            if not solve.success:
                raise RuntimeError(
                    f"local block {node}/{choice} infeasible: {solve.message}"
                )
            pseudo_result = SimpleNamespace(
                ineqlin=SimpleNamespace(
                    marginals=np.asarray(solve.x[: len(inequalities)])
                ),
                eqlin=SimpleNamespace(
                    marginals=np.asarray(solve.x[len(inequalities) :])
                ),
            )
            exact_dual = reconstruct_exact_dual(
                pseudo_result,
                inequalities,
                inequality_rhs,
                inequality_metadata,
                equalities,
                equality_rhs,
                equality_metadata,
                target,
                None,
                active_tolerance=1e-12,
            )
            if exact_dual is None:
                exact_dual = exactify_basis(
                    solve,
                    matrix,
                    target,
                    len(inequalities),
                    inequality_metadata,
                    equality_metadata,
                )
            total_inequality_rows += len(
                exact_dual["active_inequalities"]
            )
            total_equality_rows += len(exact_dual["active_equalities"])
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
                    "certificate": exact_dual,
                }
            )
        certificates[str(node)] = node_certificates
        print(
            json.dumps(
                {
                    "node": node,
                    "exact_local_blocks": len(node_certificates),
                }
            ),
            flush=True,
        )

    required_simplex_sum = (
        F(-80) - F(global_payload["global_variable_objective"])
    )
    if sum(simplex.values(), F(0)) != required_simplex_sum:
        raise RuntimeError(
            "simplex anchors do not close the global dual objective"
        )
    payload = {
        "status": "perspective_local_blocks_exact",
        "source": str(args.terminal_report),
        "global_support": str(args.global_support),
        "simplex": {
            str(node): str(value) for node, value in simplex.items()
        },
        "required_simplex_sum": str(required_simplex_sum),
        "local_block_count": sum(
            len(node_certificates)
            for node_certificates in certificates.values()
        ),
        "active_local_inequality_row_count": total_inequality_rows,
        "active_local_equality_row_count": total_equality_rows,
        "certificates": certificates,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: payload[key]
                for key in (
                    "status",
                    "simplex",
                    "required_simplex_sum",
                    "local_block_count",
                    "active_local_inequality_row_count",
                    "active_local_equality_row_count",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
