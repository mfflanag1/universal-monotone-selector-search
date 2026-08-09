#!/usr/bin/env python3
"""Certify mixed terminal topologies with a continuous perspective relaxation."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.optimize._highspy._core import (
    HighsBasisStatus,
    HighsLp,
    MatrixFormat,
    _Highs,
    kHighsInf,
)
from scipy.sparse import vstack


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
sys.path.insert(0, str(PROJECT_SRC))

from n5_facet_search import load_facets
from n5_mixed_terminal_branch_search import (
    extreme_representations,
    sparse_matrix,
)
from n5_terminal_topology_relaxation import reconstruct_exact_dual


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--float-support-output", type=Path)
    parser.add_argument("--basis-support-output", type=Path)
    parser.add_argument("--float-support-input", type=Path)
    parser.add_argument("--skip-exact-reconstruction", action="store_true")
    parser.add_argument("--variable-grand", action="store_true")
    parser.add_argument(
        "--method",
        choices=("highs-ds", "highs-ipm", "highs"),
        default="highs-ds",
        help="SciPy/HiGHS algorithm used for the floating solve",
    )
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    node_index = {node: index for index, node in enumerate(node_ids)}
    base_variable_count = len(node_ids) * coalition_count

    def column(node: int, coalition: int) -> int:
        return node_index[node] * coalition_count + coalition

    mixed = [
        terminal
        for terminal in terminals
        if terminal.get("coefficient") is None
        and len(set(terminal["scaled_divergence"])) > 2
    ]
    representations = {
        int(terminal["node"]): extreme_representations(
            [F(value) for value in terminal["scaled_divergence"]], n
        )
        for terminal in mixed
    }
    if not mixed or any(not choices for choices in representations.values()):
        raise RuntimeError("terminal report has no complete mixed representation")

    choice_offsets: dict[int, int] = {}
    copy_offsets: dict[int, int] = {}
    next_variable = base_variable_count
    for terminal in mixed:
        node = int(terminal["node"])
        choice_offsets[node] = next_variable
        next_variable += len(representations[node])
        copy_offsets[node] = next_variable
        next_variable += len(representations[node]) * coalition_count
    variable_count = next_variable

    def choice_column(node: int, choice: int) -> int:
        return choice_offsets[node] + choice

    def copy_column(node: int, choice: int, coalition: int) -> int:
        return copy_offsets[node] + choice * coalition_count + coalition

    inequalities: list[dict[int, float]] = []
    inequality_rhs: list[float] = []
    inequality_metadata: list[tuple[Any, ...]] = []
    equalities: list[dict[int, float]] = []
    equality_rhs: list[float] = []
    equality_metadata: list[tuple[Any, ...]] = []
    facets, _ = load_facets()

    for node in node_ids:
        equalities.append({column(node, 0): 1.0})
        equality_rhs.append(0.0)
        equality_metadata.append(("empty", node))
        if args.variable_grand:
            inequalities.append({column(node, grand): 1.0})
            inequality_rhs.append(1.0)
            inequality_metadata.append(("grand_upper", node))
        else:
            equalities.append({column(node, grand): 1.0})
            equality_rhs.append(1.0)
            equality_metadata.append(("grand", node))
        for facet_index, facet in enumerate(facets):
            inequalities.append(
                {
                    column(node, coalition): -float(coefficient)
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                }
            )
            inequality_rhs.append(0.0)
            inequality_metadata.append(("exact_cone", node, facet_index))
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
                inequality_metadata.append(
                    ("game_monotonicity", node, coalition, successor)
                )

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
                equality_metadata.append(
                    ("protected_invariance", source, sink, coalition)
                )
            else:
                inequalities.append(row)
                inequality_rhs.append(0.0)
                inequality_metadata.append(
                    ("directed_monotonicity", source, sink, coalition)
                )

    for terminal in mixed:
        node = int(terminal["node"])
        choices = representations[node]
        equalities.append(
            {
                choice_column(node, choice): 1.0
                for choice in range(len(choices))
            }
        )
        equality_rhs.append(1.0)
        equality_metadata.append(("perspective_simplex", node))
        for choice in range(len(choices)):
            inequalities.append({choice_column(node, choice): -1.0})
            inequality_rhs.append(0.0)
            inequality_metadata.append(
                ("perspective_choice_nonnegative", node, choice)
            )
        for coalition in range(coalition_count):
            row = {column(node, coalition): 1.0}
            for choice in range(len(choices)):
                row[copy_column(node, choice, coalition)] = -1.0
            equalities.append(row)
            equality_rhs.append(0.0)
            equality_metadata.append(
                ("perspective_link", node, coalition)
            )
        for choice in range(len(choices)):
            equalities.append({copy_column(node, choice, 0): 1.0})
            equality_rhs.append(0.0)
            equality_metadata.append(
                ("perspective_copy_empty", node, choice)
            )
            grand_row = {
                copy_column(node, choice, grand): 1.0,
                choice_column(node, choice): -1.0,
            }
            if args.variable_grand:
                inequalities.append(grand_row)
                inequality_rhs.append(0.0)
                inequality_metadata.append(
                    ("perspective_copy_grand_upper", node, choice)
                )
            else:
                equalities.append(grand_row)
                equality_rhs.append(0.0)
                equality_metadata.append(
                    ("perspective_copy_grand", node, choice)
                )
            for facet_index, facet in enumerate(facets):
                inequalities.append(
                    {
                        copy_column(node, choice, coalition): -float(
                            coefficient
                        )
                        for coalition, coefficient in enumerate(facet)
                        if coefficient
                    }
                )
                inequality_rhs.append(0.0)
                inequality_metadata.append(
                    (
                        "perspective_copy_exact_cone",
                        node,
                        choice,
                        facet_index,
                    )
                )
            for coalition in range(coalition_count):
                for player in range(n):
                    if coalition >> player & 1:
                        continue
                    successor = coalition | (1 << player)
                    inequalities.append(
                        {
                            copy_column(node, choice, coalition): 1.0,
                            copy_column(node, choice, successor): -1.0,
                        }
                    )
                    inequality_rhs.append(0.0)
                    inequality_metadata.append(
                        (
                            "perspective_copy_game_monotonicity",
                            node,
                            choice,
                            coalition,
                            successor,
                        )
                    )

    objective_exact = [F(0)] * variable_count
    objective_constant = F(0)
    mixed_nodes = set(representations)
    for terminal in terminals:
        node = int(terminal["node"])
        if node in mixed_nodes:
            continue
        raw_coefficient = terminal.get("coefficient")
        if raw_coefficient is not None:
            coefficient = F(raw_coefficient)
            coalition = int(terminal["coalition"])
            if coefficient > 0:
                objective_exact[column(node, coalition)] -= coefficient
            elif coefficient < 0:
                mass = -coefficient
                objective_exact[column(node, grand ^ coalition)] -= mass
                if args.variable_grand:
                    objective_exact[column(node, grand)] += mass
                else:
                    objective_constant -= mass
        else:
            divergence = [F(value) for value in terminal["scaled_divergence"]]
            low, high = sorted(set(divergence))
            coalition = sum(
                1 << player
                for player, value in enumerate(divergence)
                if value == high
            )
            objective_exact[column(node, coalition)] -= high - low
            objective_exact[column(node, grand)] -= low
    for node, choices in representations.items():
        for choice, (mu, coalitions) in enumerate(choices):
            if args.variable_grand:
                objective_exact[
                    copy_column(node, choice, grand)
                ] -= mu
            else:
                objective_exact[choice_column(node, choice)] -= mu
            for coalition, weight in coalitions:
                objective_exact[
                    copy_column(node, choice, coalition)
                ] -= weight

    objective = np.asarray([float(value) for value in objective_exact])
    if args.float_support_input is None:
        result = linprog(
            objective,
            A_ub=sparse_matrix(inequalities, variable_count),
            b_ub=np.asarray(inequality_rhs),
            A_eq=sparse_matrix(equalities, variable_count),
            b_eq=np.asarray(equality_rhs),
            bounds=[(None, None)] * variable_count,
            method=args.method,
            options={
                "dual_feasibility_tolerance": 1e-10,
                "primal_feasibility_tolerance": 1e-10,
            },
        )
        if not result.success:
            raise RuntimeError(result.message)
    else:
        saved_support = json.loads(args.float_support_input.read_text())
        inequality_lookup = {
            tuple(metadata): index
            for index, metadata in enumerate(inequality_metadata)
        }
        equality_lookup = {
            tuple(metadata): index
            for index, metadata in enumerate(equality_metadata)
        }
        inequality_marginals = np.zeros(len(inequalities))
        equality_marginals = np.zeros(len(equalities))
        for active in saved_support["inequalities"]:
            inequality_marginals[
                inequality_lookup[tuple(active["row"])]
            ] = float(active["weight"])
        for active in saved_support["equalities"]:
            equality_marginals[
                equality_lookup[tuple(active["row"])]
            ] = float(active["weight"])
        result = SimpleNamespace(
            fun=float(objective_constant),
            ineqlin=SimpleNamespace(marginals=inequality_marginals),
            eqlin=SimpleNamespace(marginals=equality_marginals),
        )
    exact_dual = (
        None
        if args.skip_exact_reconstruction
        else reconstruct_exact_dual(
            result,
            inequalities,
            inequality_rhs,
            inequality_metadata,
            equalities,
            equality_rhs,
            equality_metadata,
            objective_exact,
            objective_constant,
            active_tolerance=1e-12,
        )
    )
    maximum = -float(result.fun) + float(objective_constant)
    nonzero_inequality_marginals = [
        float(value)
        for value in result.ineqlin.marginals
        if abs(float(value)) > 1e-12
    ]
    nonzero_equality_marginals = [
        float(value)
        for value in result.eqlin.marginals
        if abs(float(value)) > 1e-12
    ]
    payload: dict[str, Any] = {
        "status": (
            "mixed_terminal_perspective_relaxation_exact"
            if exact_dual is not None
            else "mixed_terminal_perspective_relaxation_float"
        ),
        "source": str(args.terminal_report),
        "variable_grand": args.variable_grand,
        "solver_method": args.method,
        "terminal_count": len(terminals),
        "protected_terminal_pair_count": len(protected_players),
        "mixed_terminal_count": len(mixed),
        "representation_counts": {
            str(node): len(choices)
            for node, choices in representations.items()
        },
        "variable_count": variable_count,
        "inequality_count": len(inequalities),
        "equality_count": len(equalities),
        "maximum_lower_expectation_objective_float": maximum,
        "nonzero_inequality_dual_count": len(
            nonzero_inequality_marginals
        ),
        "nonzero_equality_dual_count": len(
            nonzero_equality_marginals
        ),
        "dual_reconstruction_max_denominator_at_1e6": max(
            (
                F(value).limit_denominator(1_000_000).denominator
                for value in (
                    nonzero_inequality_marginals
                    + nonzero_equality_marginals
                )
            ),
            default=1,
        ),
        "safe_for_every_terminal_realization": (
            maximum <= 1e-9 and exact_dual is not None
        ),
        "exact_dual_verified": exact_dual is not None,
        "exact_dual_certificate": exact_dual,
        "interpretation": (
            "The perspective variables continuously convexify every extreme "
            "core-dual branch at each genuinely mixed terminal. Every original "
            "terminal realization is feasible in this relaxation, so an exact "
            "nonpositive optimum is a topology-level no-go proof."
        ),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    if args.float_support_output is not None:
        support_payload = {
            "inequalities": [
                {
                    "row": list(inequality_metadata[index]),
                    "weight": format(float(value), ".17g"),
                }
                for index, value in enumerate(result.ineqlin.marginals)
                if abs(float(value)) > 1e-12
            ],
            "equalities": [
                {
                    "row": list(equality_metadata[index]),
                    "weight": format(float(value), ".17g"),
                }
                for index, value in enumerate(result.eqlin.marginals)
                if abs(float(value)) > 1e-12
            ],
        }
        args.float_support_output.write_text(
            json.dumps(support_payload, indent=2) + "\n"
        )
    if args.basis_support_output is not None:
        inequality_matrix = sparse_matrix(inequalities, variable_count)
        equality_matrix = sparse_matrix(equalities, variable_count)
        matrix = vstack((inequality_matrix, equality_matrix)).tocsc()
        lp = HighsLp()
        lp.num_col_ = variable_count
        lp.num_row_ = matrix.shape[0]
        lp.col_cost_ = objective
        lp.col_lower_ = np.full(variable_count, -kHighsInf)
        lp.col_upper_ = np.full(variable_count, kHighsInf)
        lp.row_lower_ = np.concatenate(
            (
                np.full(len(inequalities), -kHighsInf),
                np.asarray(equality_rhs),
            )
        )
        lp.row_upper_ = np.concatenate(
            (
                np.asarray(inequality_rhs),
                np.asarray(equality_rhs),
            )
        )
        lp.a_matrix_.format_ = MatrixFormat.kColwise
        lp.a_matrix_.num_col_ = variable_count
        lp.a_matrix_.num_row_ = matrix.shape[0]
        lp.a_matrix_.start_ = matrix.indptr.astype(np.int32)
        lp.a_matrix_.index_ = matrix.indices.astype(np.int32)
        lp.a_matrix_.value_ = matrix.data
        highs = _Highs()
        highs.setOptionValue("output_flag", False)
        highs.setOptionValue("solver", "simplex")
        highs.setOptionValue("dual_feasibility_tolerance", 1e-10)
        highs.setOptionValue("primal_feasibility_tolerance", 1e-10)
        if highs.passModel(lp).name != "kOk":
            raise RuntimeError("HiGHS rejected the perspective model")
        if highs.run().name != "kOk":
            raise RuntimeError("HiGHS failed on the perspective model")
        basis = highs.getBasis()
        solution = highs.getSolution()
        nonbasic_rows = [
            index
            for index, status in enumerate(basis.row_status)
            if status != HighsBasisStatus.kBasic
        ]
        if len(nonbasic_rows) != variable_count:
            raise RuntimeError(
                f"basis has {len(nonbasic_rows)} nonbasic rows for "
                f"{variable_count} free columns"
            )
        basis_support = {
            "status": "perspective_simplex_basis_support",
            "nonbasic_row_count": len(nonbasic_rows),
            "variable_count": variable_count,
            "inequalities": [
                {
                    "row": list(inequality_metadata[index]),
                    "weight": format(float(solution.row_dual[index]), ".17g"),
                }
                for index in nonbasic_rows
                if index < len(inequalities)
            ],
            "equalities": [
                {
                    "row": list(
                        equality_metadata[index - len(inequalities)]
                    ),
                    "weight": format(float(solution.row_dual[index]), ".17g"),
                }
                for index in nonbasic_rows
                if index >= len(inequalities)
            ],
        }
        args.basis_support_output.write_text(
            json.dumps(basis_support, indent=2) + "\n"
        )
    print(
        json.dumps(
            {
                key: payload[key]
                for key in (
                    "status",
                    "variable_count",
                    "inequality_count",
                    "equality_count",
                    "maximum_lower_expectation_objective_float",
                    "safe_for_every_terminal_realization",
                    "exact_dual_verified",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
