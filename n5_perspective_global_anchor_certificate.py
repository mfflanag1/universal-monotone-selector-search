#!/usr/bin/env python3
"""Exactify the global terminal block for a rational perspective anchor."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack

from n5_facet_search import load_facets
from n5_mixed_terminal_branch_search import sparse_matrix
from n5_terminal_topology_relaxation import reconstruct_exact_dual


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("anchor_search", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--float-support-output", type=Path)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    anchor_payload = json.loads(args.anchor_search.read_text())
    anchor = anchor_payload["best_complete_anchor"]
    if anchor is None:
        raise RuntimeError("anchor search has no complete rational anchor")
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    node_index = {node: index for index, node in enumerate(node_ids)}
    variable_count = len(node_ids) * coalition_count

    def column(node: int, coalition: int) -> int:
        return node_index[node] * coalition_count + coalition

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

    links = {
        int(node): [F(value) for value in values]
        for node, values in anchor["links"].items()
    }
    mixed_nodes = set(links)
    objective_exact = [F(0)] * variable_count
    objective_constant = F(0)
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
    for node, values in links.items():
        for coalition, value in enumerate(values):
            objective_exact[column(node, coalition)] -= value

    dual_matrix = sparse_matrix(
        inequalities + equalities, variable_count
    ).transpose().tocsr()
    dual_rhs = np.asarray(inequality_rhs + equality_rhs)
    simplex = {
        int(node): F(value) for node, value in anchor["simplex"].items()
    }
    target_dual_objective = objective_constant - sum(
        simplex.values(), F(0)
    )
    augmented_matrix = vstack(
        (
            dual_matrix,
            csr_matrix(dual_rhs.reshape(1, -1)),
        ),
        format="csr",
    )
    solve = linprog(
        np.zeros(len(dual_rhs)),
        A_eq=augmented_matrix,
        b_eq=np.asarray(
            [float(value) for value in objective_exact]
            + [float(target_dual_objective)]
        ),
        bounds=(
            [(None, 0.0)] * len(inequalities)
            + [(None, None)] * len(equalities)
        ),
        method="highs-ds",
        options={
            "dual_feasibility_tolerance": 1e-10,
            "primal_feasibility_tolerance": 1e-10,
        },
    )
    if not solve.success:
        raise RuntimeError(solve.message)
    result = SimpleNamespace(
        ineqlin=SimpleNamespace(
            marginals=np.asarray(solve.x[: len(inequalities)])
        ),
        eqlin=SimpleNamespace(
            marginals=np.asarray(solve.x[len(inequalities) :])
        ),
    )
    exact_dual = reconstruct_exact_dual(
        result,
        inequalities,
        inequality_rhs,
        inequality_metadata,
        equalities,
        equality_rhs,
        equality_metadata,
        objective_exact,
        target_dual_objective,
        active_tolerance=1e-12,
    )
    if exact_dual is None:
        if args.float_support_output is not None:
            args.float_support_output.write_text(
                json.dumps(
                    {
                        "inequalities": [
                            {
                                "row": list(inequality_metadata[index]),
                                "weight": format(float(value), ".17g"),
                            }
                            for index, value in enumerate(
                                result.ineqlin.marginals
                            )
                            if abs(float(value)) > 1e-12
                        ],
                        "equalities": [
                            {
                                "row": list(equality_metadata[index]),
                                "weight": format(float(value), ".17g"),
                            }
                            for index, value in enumerate(
                                result.eqlin.marginals
                            )
                            if abs(float(value)) > 1e-12
                        ],
                    },
                    indent=2,
                )
                + "\n"
            )
        raise RuntimeError("failed to reconstruct the global anchor dual")
    base_variable_objective = F(exact_dual["variable_objective"])
    adjusted_simplex = dict(simplex)
    adjusted_node = max(adjusted_simplex)
    adjusted_simplex[adjusted_node] = (
        objective_constant
        - base_variable_objective
        - sum(
            (
                value
                for node, value in adjusted_simplex.items()
                if node != adjusted_node
            ),
            F(0),
        )
    )
    payload = {
        "status": "perspective_global_anchor_exact",
        "source": str(args.terminal_report),
        "anchor_search": str(args.anchor_search),
        "objective_constant": str(objective_constant),
        "base_variable_objective": str(base_variable_objective),
        "original_simplex": {
            str(node): str(value) for node, value in simplex.items()
        },
        "adjusted_simplex": {
            str(node): str(value)
            for node, value in adjusted_simplex.items()
        },
        "simplex_adjustment": str(
            adjusted_simplex[adjusted_node] - simplex[adjusted_node]
        ),
        "links": {
            str(node): [str(value) for value in values]
            for node, values in links.items()
        },
        "exact_base_dual_certificate": exact_dual,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: payload[key]
                for key in (
                    "status",
                    "objective_constant",
                    "base_variable_objective",
                    "adjusted_simplex",
                    "simplex_adjustment",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
