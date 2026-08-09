#!/usr/bin/env python3
"""Compress a perspective dual into exact terminal-envelope certificates."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
import z3
from scipy.linalg import null_space, qr
from scipy.optimize import linprog
from scipy.optimize._highspy._core import (
    HighsBasisStatus,
    HighsLp,
    MatrixFormat,
    _Highs,
    kHighsInf,
)
from scipy.sparse import coo_matrix, vstack


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
# Prefer the auditable sources in this worktree while using the research
# project only for modules (notably the facet catalogue) absent locally.
sys.path.append(str(PROJECT_SRC))

from n5_middle_split_identity_certificate import (
    coefficient_vector,
    exact_generators,
    exact_solution_on_support,
    find_certificate,
)
from n5_mixed_divergence_premium import choquet_representation
from n5_mixed_terminal_branch_search import extreme_representations


F = Fraction


class TopologyReconstructionError(RuntimeError):
    """A feasible float topology LP whose current basis did not exactify."""


def exact_topology_projection(
    envelopes: dict[int, tuple[F, ...]],
    terminal_divergences: dict[int, tuple[F, ...]],
    node_ids: list[int],
    n: int,
    generators: list[tuple[F, ...]],
    epsilon: F,
    envelope_cuts: list[tuple[int, tuple[F, ...], F]],
    branch_memberships: list[tuple[int, int, tuple[F, ...]]],
    local_cone_generators: list[tuple[F, ...]],
    topology_objective_salt: int = 0,
) -> tuple[dict[int, tuple[F, ...]], list[tuple[int, F]]]:
    """Project nearby envelopes onto the topology dual cone exactly."""
    coalition_count = 1 << n
    topology_count = len(generators)
    local_count = len(local_cone_generators)
    variable_count = topology_count + len(branch_memberships) * local_count
    node_index = {node: index for index, node in enumerate(node_ids)}

    topology_columns_by_coordinate = [
        [
            (index, generator[coordinate])
            for index, generator in enumerate(generators)
            if generator[coordinate]
        ]
        for coordinate in range(len(node_ids) * coalition_count)
    ]
    local_columns_by_coordinate = [
        [
            (index, generator[coordinate])
            for index, generator in enumerate(local_cone_generators)
            if generator[coordinate]
        ]
        for coordinate in range(coalition_count)
    ]
    equality_rows: list[dict[int, F]] = []
    equality_rhs: list[F] = []
    inequality_rows: list[dict[int, F]] = []
    inequality_rhs: list[F] = []

    for node_offset, node in enumerate(node_ids):
        for coalition in range(coalition_count):
            coordinate = node_offset * coalition_count + coalition
            expected = envelopes[node][coalition]
            # envelope = - topology_generator_matrix * topology_weights
            inequality_rows.append(
                {
                    index: -value
                    for index, value in topology_columns_by_coordinate[
                        coordinate
                    ]
                }
            )
            inequality_rhs.append(expected + epsilon)
            inequality_rows.append(
                dict(topology_columns_by_coordinate[coordinate])
            )
            inequality_rhs.append(-(expected - epsilon))
        for player, divergence_value in enumerate(
            terminal_divergences[node]
        ):
            row: dict[int, F] = defaultdict(F)
            for coalition in range(coalition_count):
                if not coalition >> player & 1:
                    continue
                coordinate = node_offset * coalition_count + coalition
                for index, value in topology_columns_by_coordinate[
                    coordinate
                ]:
                    row[index] -= value
            equality_rows.append(dict(row))
            equality_rhs.append(divergence_value)
    topology_equality_count = len(equality_rows)
    for node, game, branch_value in envelope_cuts:
        row: dict[int, F] = defaultdict(F)
        offset = node_index[node] * coalition_count
        for coalition, game_value in enumerate(game):
            if not game_value:
                continue
            for index, value in topology_columns_by_coordinate[
                offset + coalition
            ]:
                row[index] += value * game_value
        inequality_rows.append(dict(row))
        inequality_rhs.append(-branch_value)
    for membership_index, (node, _, branch) in enumerate(
        branch_memberships
    ):
        local_offset = topology_count + membership_index * local_count
        node_offset = node_index[node] * coalition_count
        for coalition in range(coalition_count):
            row: dict[int, F] = defaultdict(F)
            for index, value in topology_columns_by_coordinate[
                node_offset + coalition
            ]:
                row[index] -= value
            for index, value in local_columns_by_coordinate[coalition]:
                row[local_offset + index] -= value
            equality_rows.append(dict(row))
            equality_rhs.append(branch[coalition])

    def sparse_rows(
        rows: list[dict[int, F]],
        column_count: int = variable_count,
        row_scales: list[F] | None = None,
    ):
        row_indices = []
        column_indices = []
        values = []
        for row_index, row in enumerate(rows):
            scale = F(1) if row_scales is None else row_scales[row_index]
            for column, value in row.items():
                if value:
                    row_indices.append(row_index)
                    column_indices.append(column)
                    values.append(float(value / scale))
        return coo_matrix(
            (values, (row_indices, column_indices)),
            shape=(len(rows), column_count),
        ).tocsc()

    objective = np.asarray(
        [1.0 + (index % 97) * 1e-7 for index in range(variable_count)]
    )
    inequality_scales = [
        max(
            [F(1), abs(rhs)]
            + [abs(value) for value in row.values()]
        )
        for row, rhs in zip(
            inequality_rows, inequality_rhs, strict=True
        )
    ]
    equality_scales = [
        max(
            [F(1), abs(rhs)]
            + [abs(value) for value in row.values()]
        )
        for row, rhs in zip(equality_rows, equality_rhs, strict=True)
    ]
    float_inequality_matrix = sparse_rows(
        inequality_rows, row_scales=inequality_scales
    )
    float_inequality_rhs = np.asarray(
        [
            float(value / scale)
            for value, scale in zip(
                inequality_rhs, inequality_scales, strict=True
            )
        ]
    )
    float_equality_matrix = sparse_rows(
        equality_rows, row_scales=equality_scales
    )
    float_equality_rhs = np.asarray(
        [
            float(value / scale)
            for value, scale in zip(
                equality_rhs, equality_scales, strict=True
            )
        ]
    )
    floated = None
    floating_failures = []
    for method, options in (
        (
            "highs-ds",
            {
                "dual_feasibility_tolerance": 1e-10,
                "primal_feasibility_tolerance": 1e-10,
            },
        ),
        (
            "highs-ipm",
            {
                "dual_feasibility_tolerance": 1e-9,
                "primal_feasibility_tolerance": 1e-9,
                "ipm_optimality_tolerance": 1e-10,
            },
        ),
        ("highs", {}),
    ):
        candidate = linprog(
            objective,
            A_ub=float_inequality_matrix,
            b_ub=float_inequality_rhs,
            A_eq=float_equality_matrix,
            b_eq=float_equality_rhs,
            bounds=(0.0, None),
            method=method,
            options=options,
        )
        if candidate.success:
            floated = candidate
            if floating_failures:
                print(
                    json.dumps(
                        {
                            "topology_projection_float_recovery_method": (
                                method
                            ),
                            "prior_failure_count": len(floating_failures),
                        }
                    ),
                    flush=True,
                )
            break
        floating_failures.append((method, candidate.message))
    if floated is None:
        raise RuntimeError(
            "floating topology projection failed with every algorithm: "
            + "; ".join(
                f"{method}: {message}"
                for method, message in floating_failures
            )
        )

    # Keep the topology-only exactification close to the topology component of
    # the membership-aware solve.  Without these narrow internal guide boxes,
    # re-optimizing only the topology weights can move to a distant point in
    # the much wider rounding box and discard all of the local-membership
    # information that the first LP supplied.  These rows are not assumptions
    # in the archived proof: generate() still rebuilds every local membership
    # exactly and verify() replays the resulting identities independently.
    topology_anchor_radius = F(1, 100_000_000)
    floated_topology = floated.x[:topology_count]
    for node_offset, _ in enumerate(node_ids):
        for coalition in range(coalition_count):
            coordinate = node_offset * coalition_count + coalition
            anchor = F(
                -sum(
                    float(value) * float(floated_topology[index])
                    for index, value in topology_columns_by_coordinate[
                        coordinate
                    ]
                )
            ).limit_denominator(1_000_000_000_000)
            inequality_rows.append(
                {
                    index: -value
                    for index, value in topology_columns_by_coordinate[
                        coordinate
                    ]
                }
            )
            inequality_rhs.append(anchor + topology_anchor_radius)
            inequality_rows.append(
                dict(topology_columns_by_coordinate[coordinate])
            )
            inequality_rhs.append(-anchor + topology_anchor_radius)

    # The local variables are feasibility witnesses used to steer the float
    # projection.  They are not part of the global certificate: local cone
    # membership is reconstructed independently (and exactly) below.  Exactify
    # only the topology weights against the topology equalities and all cuts;
    # any lost local membership produces an exact separating cut in generate().
    exact_equality_rows = equality_rows[:topology_equality_count]
    exact_equality_rhs = equality_rhs[:topology_equality_count]
    topology_objective = objective[:topology_count]
    if topology_objective_salt:
        topology_indices = np.arange(topology_count, dtype=np.int64)
        topology_objective = 1.0 + (
            (
                topology_indices * (2 * topology_objective_salt + 1)
                + topology_objective_salt * topology_objective_salt
            )
            % 1009
        ) * 1e-6
    topology_inequality_scales = [
        max(
            [F(1), abs(rhs)]
            + [abs(value) for value in row.values()]
        )
        for row, rhs in zip(
            inequality_rows, inequality_rhs, strict=True
        )
    ]
    topology_inequality_matrix = sparse_rows(
        inequality_rows,
        topology_count,
        topology_inequality_scales,
    )
    topology_equality_matrix = sparse_rows(
        exact_equality_rows, topology_count
    )
    topology_floated = None
    topology_replay_margin = F(0)
    working_inequality_rhs = inequality_rhs
    maximum_margin = min(epsilon / 1_000, F(1, 1_000_000_000))
    for candidate_margin in (
        maximum_margin,
        maximum_margin / 10,
        maximum_margin / 100,
        F(0),
    ):
        candidate_rhs = [
            value - candidate_margin for value in inequality_rhs
        ]
        candidate_solution = linprog(
            topology_objective,
            A_ub=topology_inequality_matrix,
            b_ub=np.asarray(
                [
                    float(value / scale)
                    for value, scale in zip(
                        candidate_rhs,
                        topology_inequality_scales,
                        strict=True,
                    )
                ]
            ),
            A_eq=topology_equality_matrix,
            b_eq=np.asarray([float(value) for value in exact_equality_rhs]),
            bounds=(0.0, None),
            method="highs-ds",
            options={
                "dual_feasibility_tolerance": 1e-10,
                "primal_feasibility_tolerance": 1e-10,
            },
        )
        if candidate_solution.success:
            topology_floated = candidate_solution
            topology_replay_margin = candidate_margin
            working_inequality_rhs = candidate_rhs
            break
    if topology_floated is None:
        raise RuntimeError(
            "floating topology-only projection failed at every safety margin"
        )
    topology_matrix = vstack(
        [topology_inequality_matrix, topology_equality_matrix],
        format="csc",
    )
    topology_lp = HighsLp()
    topology_lp.num_col_ = topology_count
    topology_lp.num_row_ = len(inequality_rows) + len(exact_equality_rows)
    topology_lp.col_cost_ = topology_objective
    topology_lp.col_lower_ = np.zeros(topology_count)
    topology_lp.col_upper_ = np.full(topology_count, kHighsInf)
    topology_lp.row_lower_ = np.concatenate(
        [
            np.full(len(inequality_rows), -kHighsInf),
            np.asarray([float(value) for value in exact_equality_rhs]),
        ]
    )
    topology_lp.row_upper_ = np.concatenate(
        [
            np.asarray(
                [
                    float(value / scale)
                    for value, scale in zip(
                        working_inequality_rhs,
                        topology_inequality_scales,
                        strict=True,
                    )
                ]
            ),
            np.asarray([float(value) for value in exact_equality_rhs]),
        ]
    )
    topology_lp.a_matrix_.format_ = MatrixFormat.kColwise
    topology_lp.a_matrix_.num_col_ = topology_count
    topology_lp.a_matrix_.num_row_ = topology_lp.num_row_
    topology_lp.a_matrix_.start_ = topology_matrix.indptr.astype(np.int32)
    topology_lp.a_matrix_.index_ = topology_matrix.indices.astype(np.int32)
    topology_lp.a_matrix_.value_ = topology_matrix.data
    topology_highs = _Highs()
    topology_highs.setOptionValue("output_flag", False)
    topology_highs.setOptionValue("solver", "simplex")
    topology_highs.setOptionValue("dual_feasibility_tolerance", 1e-10)
    topology_highs.setOptionValue("primal_feasibility_tolerance", 1e-10)
    topology_pass_status = topology_highs.passModel(topology_lp)
    topology_run_status = topology_highs.run()
    if (
        topology_pass_status.name not in {"kOk", "kWarning"}
        or topology_run_status.name != "kOk"
        or topology_highs.getModelStatus().name != "kOptimal"
    ):
        raise RuntimeError(
            "direct topology basis solve failed: "
            f"pass={topology_pass_status.name}, "
            f"run={topology_run_status.name}, "
            f"model={topology_highs.getModelStatus().name}"
        )
    topology_basis = topology_highs.getBasis()
    topology_solution = topology_highs.getSolution()
    topology_x = np.asarray(topology_solution.col_value, dtype=float)
    basic_topology_variables = {
        index
        for index, status in enumerate(topology_basis.col_status)
        if status == HighsBasisStatus.kBasic
    }
    inequality_slacks = [
        (
            float(rhs)
            - sum(
                float(coefficient) * float(topology_x[column])
                for column, coefficient in row.items()
            )
        )
        / float(scale)
        for row, rhs, scale in zip(
            inequality_rows,
            working_inequality_rhs,
            topology_inequality_scales,
            strict=True,
        )
    ]
    dual_marginals = [
        float(value)
        for value in topology_solution.row_dual[: len(inequality_rows)]
    ]
    basis_active_inequalities = tuple(
        index
        for index, status in enumerate(
            topology_basis.row_status[: len(inequality_rows)]
        )
        if status != HighsBasisStatus.kBasic
    )
    base_active_candidates: list[tuple[str, tuple[int, ...]]] = [
        ("simplex_row_basis", basis_active_inequalities)
    ]
    base_active_candidates.extend(
        (
            f"dual_marginal_{tolerance:g}",
            tuple(
                index
                for index, marginal in enumerate(dual_marginals)
                if abs(marginal) > tolerance
            ),
        )
        for tolerance in (1e-7, 1e-9, 1e-11, 1e-13, 0.0)
    )
    base_active_candidates.extend(
        (
            f"slack_{tolerance:g}",
            tuple(
                index
                for index, slack in enumerate(inequality_slacks)
                if slack <= tolerance
            ),
        )
        for tolerance in (1e-9, 1e-7, 1e-6)
    )
    base_active_candidates.append(("none", ()))
    tried_systems: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    for threshold in (1e-7, 1e-9, 1e-11, 1e-13, 0.0):
        variable_support = tuple(
            index
            for index, value in enumerate(topology_x)
            if index in basic_topology_variables or value > threshold
        )
        near_active = tuple(
            index
            for index, slack in enumerate(inequality_slacks)
            if slack <= 1e-9
        )
        dense_equality_matrix = np.asarray(
            [
                [float(row.get(variable, F(0))) for variable in variable_support]
                for row in exact_equality_rows
            ],
            dtype=float,
        )
        equality_kernel = null_space(dense_equality_matrix, rcond=1e-11)
        dense_active_matrix = np.asarray(
            [
                [
                    float(inequality_rows[index].get(variable, F(0)))
                    for variable in variable_support
                ]
                for index in near_active
            ],
            dtype=float,
        )
        projected_active_matrix = dense_active_matrix @ equality_kernel
        row_norms = np.linalg.norm(projected_active_matrix, axis=1)
        scaled_rank_matrix = projected_active_matrix.copy()
        nonzero_rows = row_norms > 0
        scaled_rank_matrix[nonzero_rows] /= row_norms[nonzero_rows, None]
        _, triangular, pivots = qr(
            scaled_rank_matrix.T, mode="economic", pivoting=True
        )
        diagonal = np.abs(np.diag(triangular))
        numerical_rank = int(
            np.sum(
                diagonal
                > (
                    max(scaled_rank_matrix.shape)
                    * np.finfo(float).eps
                    * (float(diagonal[0]) if len(diagonal) else 1.0)
                )
            )
        )
        selected_rank_rows = set(
            int(index) for index in pivots[:numerical_rank]
        )
        rank_active_inequalities = tuple(
            near_active[index]
            for index in sorted(selected_rank_rows)
        )
        equality_rank = len(variable_support) - equality_kernel.shape[1]
        active_candidates = [
            (
                f"conditional_rank_basis_{equality_rank + numerical_rank}",
                rank_active_inequalities,
            ),
            *base_active_candidates,
        ]
        for active_selector, active_inequalities in active_candidates:
            system_key = (variable_support, active_inequalities)
            if system_key in tried_systems:
                continue
            tried_systems.add(system_key)

            selected_rows = exact_equality_rows + [
                inequality_rows[index]
                for index in active_inequalities
            ]
            target = tuple(exact_equality_rhs) + tuple(
                working_inequality_rhs[index]
                for index in active_inequalities
            )
            columns = [
                tuple(row.get(variable, F(0)) for row in selected_rows)
                for variable in variable_support
            ]
            print(
                json.dumps(
                    {
                        "topology_projection_reconstruction_attempt": True,
                        "variable_support_count": len(variable_support),
                        "equality_count": len(exact_equality_rows),
                        "active_inequality_count": len(
                            active_inequalities
                        ),
                        "support_threshold": threshold,
                        "active_inequality_selector": active_selector,
                    }
                ),
                flush=True,
            )
            compressed_support = exact_solution_on_support(
                target,
                columns,
                list(range(len(columns))),
                [
                    float(topology_x[index])
                    for index in variable_support
                ],
                (1_000_000_000_000,),
            )
            if compressed_support is None:
                print(
                    json.dumps(
                        {
                            "topology_projection_reconstruction_result": (
                                "no_nonnegative_exact_solution"
                            ),
                            "active_inequality_selector": active_selector,
                        }
                    ),
                    flush=True,
                )
                continue

            exact_x = [F(0)] * topology_count
            for compressed_index, weight in compressed_support:
                exact_x[variable_support[compressed_index]] = weight
            if any(weight < 0 for weight in exact_x):
                continue
            equality_failures = sum(
                sum(
                    (
                        coefficient * exact_x[column]
                        for column, coefficient in row.items()
                    ),
                    F(0),
                )
                != rhs
                for row, rhs in zip(
                    exact_equality_rows, exact_equality_rhs, strict=True
                )
            )
            if equality_failures:
                print(
                    json.dumps(
                        {
                            "topology_projection_reconstruction_result": (
                                "exact_equality_failure"
                            ),
                            "failure_count": equality_failures,
                            "active_inequality_selector": active_selector,
                        }
                    ),
                    flush=True,
                )
                continue
            inequality_violations = [
                sum(
                    (
                        coefficient * exact_x[column]
                        for column, coefficient in row.items()
                    ),
                    F(0),
                )
                - rhs
                for row, rhs in zip(
                    inequality_rows, inequality_rhs, strict=True
                )
            ]
            positive_violations = [
                value for value in inequality_violations if value > 0
            ]
            # These inequalities guide the projection but are not premises of
            # the archived proof.  A tiny simplex-basis perturbation is safe
            # here because generate() subsequently reconstructs every local
            # cone identity exactly and verify() replays those identities plus
            # the exact global topology certificate with no tolerance.
            projection_replay_tolerance = min(
                epsilon / 1_000, F(1, 1_000_000_000)
            )
            if (
                positive_violations
                and max(positive_violations) > projection_replay_tolerance
            ):
                print(
                    json.dumps(
                        {
                            "topology_projection_reconstruction_result": (
                                "exact_inequality_failure"
                            ),
                            "failure_count": len(positive_violations),
                            "maximum_violation": str(max(positive_violations)),
                            "active_inequality_selector": active_selector,
                        }
                    ),
                    flush=True,
                )
                continue
            if positive_violations:
                print(
                    json.dumps(
                        {
                            "topology_projection_reconstruction_result": (
                                "accepted_internal_guidance_perturbation"
                            ),
                            "failure_count": len(positive_violations),
                            "maximum_violation": str(max(positive_violations)),
                            "guidance_tolerance": str(
                                projection_replay_tolerance
                            ),
                            "active_inequality_selector": active_selector,
                        }
                    ),
                    flush=True,
                )

            support = [
                (index, exact_x[index])
                for index in range(topology_count)
                if exact_x[index]
            ]
            projected: dict[int, tuple[F, ...]] = {}
            for node_offset, node in enumerate(node_ids):
                projected[node] = tuple(
                    -sum(
                        (
                            weight
                            * generators[index][
                                node_offset * coalition_count + coalition
                            ]
                            for index, weight in support
                        ),
                        F(0),
                    )
                    for coalition in range(coalition_count)
                )
            print(
                json.dumps(
                    {
                        "topology_projection_float_variable_count": (
                            variable_count
                        ),
                        "topology_projection_float_topology_variable_count": (
                            topology_count
                        ),
                        "topology_projection_exact_topology_support": len(
                            support
                        ),
                        "topology_projection_exact_local_support": sum(
                            bool(floated.x[index])
                            for index in range(topology_count, variable_count)
                        ),
                        "topology_projection_support_threshold": threshold,
                        "topology_projection_active_inequality_selector": (
                            active_selector
                        ),
                        "topology_projection_active_inequality_count": len(
                            active_inequalities
                        ),
                        "topology_projection_replay_margin": str(
                            topology_replay_margin
                        ),
                        "topology_projection_anchor_radius": str(
                            topology_anchor_radius
                        ),
                        "topology_projection_objective_salt": (
                            topology_objective_salt
                        ),
                    }
                ),
                flush=True,
            )
            return projected, support
    raise TopologyReconstructionError(
        "floating topology basis did not admit an exact rational projection"
    )


def exact_separating_games(
    envelope: tuple[F, ...],
    branch: tuple[F, ...],
    generators: list[tuple[F, ...]],
    n: int,
    maximum_game_count: int = 32,
) -> list[tuple[tuple[F, ...], F]]:
    """Return exact normalized games separating a failed envelope branch."""
    coalition_count = 1 << n
    grand = coalition_count - 1
    objective = np.asarray(
        [float(branch_value - envelope_value) for branch_value, envelope_value in zip(branch, envelope, strict=True)]
    )
    inequality_matrix = -np.asarray(
        [[float(value) for value in generator] for generator in generators]
    )
    grand_equality = np.zeros((1, coalition_count))
    grand_equality[0, grand] = 1.0
    solved = linprog(
        -objective,
        A_ub=inequality_matrix,
        b_ub=np.zeros(len(generators)),
        A_eq=grand_equality,
        b_eq=np.asarray([1.0]),
        bounds=[(None, None)] * coalition_count,
        method="highs-ds",
        options={
            "dual_feasibility_tolerance": 1e-10,
            "primal_feasibility_tolerance": 1e-10,
        },
    )
    if not solved.success:
        return []

    def rational_separator(values: np.ndarray):
        for denominator in (
            1_000,
            1_000_000,
            1_000_000_000,
            1_000_000_000_000,
        ):
            game = tuple(
                F(float(value)).limit_denominator(denominator)
                for value in values
            )
            if game[grand] != 1:
                continue
            if any(
                sum(
                    (
                        coefficient * value
                        for coefficient, value in zip(
                            generator, game, strict=True
                        )
                    ),
                    F(0),
                )
                < 0
                for generator in generators
            ):
                continue
            exact_branch_value = sum(
                (
                    coefficient * value
                    for coefficient, value in zip(
                        branch, game, strict=True
                    )
                ),
                F(0),
            )
            envelope_value = sum(
                (
                    coefficient * value
                    for coefficient, value in zip(
                        envelope, game, strict=True
                    )
                ),
                F(0),
            )
            if exact_branch_value > envelope_value:
                return game, exact_branch_value
        return None

    separators = []
    seen_games = set()
    primary = rational_separator(solved.x)
    if primary is not None:
        separators.append(primary)
        seen_games.add(primary[0])
    maximum_violation = -float(solved.fun)
    if maximum_violation > 1e-10:
        violating_matrix = np.vstack([inequality_matrix, -objective])
        violating_rhs = np.concatenate(
            [
                np.zeros(len(generators)),
                np.asarray([-maximum_violation / 4]),
            ]
        )
        for salt in range(1, maximum_game_count):
            direction = np.asarray(
                [
                    (
                        ((coalition + 1) * (2 * salt + 1) + salt * salt)
                        % 101
                    )
                    - 50
                    for coalition in range(coalition_count)
                ],
                dtype=float,
            )
            direction[grand] = 0.0
            diversified = linprog(
                -direction,
                A_ub=violating_matrix,
                b_ub=violating_rhs,
                A_eq=grand_equality,
                b_eq=np.asarray([1.0]),
                bounds=[(None, None)] * coalition_count,
                method="highs-ds",
                options={
                    "dual_feasibility_tolerance": 1e-10,
                    "primal_feasibility_tolerance": 1e-10,
                },
            )
            if not diversified.success:
                continue
            separator = rational_separator(diversified.x)
            if separator is None or separator[0] in seen_games:
                continue
            separators.append(separator)
            seen_games.add(separator[0])
    if not separators:
        # A rational topology perturbation can place the target only
        # infinitesimally outside a boundary face.  HiGHS may correctly reject
        # cone membership while its normalized separator objective is too
        # small to rationalize reliably.  Decide that residual case directly
        # in exact linear arithmetic.
        game_variables = [
            z3.Real(f"separator_game_{coalition}")
            for coalition in range(coalition_count)
        ]

        def z3_rational(value: F):
            return z3.Q(value.numerator, value.denominator)

        exact_solver = z3.Solver()
        exact_solver.set(timeout=300_000)
        exact_solver.add(game_variables[grand] == 1)
        for generator in generators:
            terms = [
                game_variables[coalition] * z3_rational(coefficient)
                for coalition, coefficient in enumerate(generator)
                if coefficient
            ]
            exact_solver.add(
                (z3.Sum(terms) if terms else z3.RealVal(0)) >= 0
            )
        exact_gap_terms = [
            game_variables[coalition]
            * z3_rational(branch_value - envelope_value)
            for coalition, (branch_value, envelope_value) in enumerate(
                zip(branch, envelope, strict=True)
            )
            if branch_value != envelope_value
        ]
        exact_solver.add(
            (
                z3.Sum(exact_gap_terms)
                if exact_gap_terms
                else z3.RealVal(0)
            )
            > 0
        )
        if exact_solver.check() == z3.sat:
            model = exact_solver.model()
            game = []
            for variable in game_variables:
                value = model.eval(variable, model_completion=True)
                if not z3.is_rational_value(value):
                    raise RuntimeError(
                        "non-rational value in exact separator model"
                    )
                game.append(
                    F(
                        value.numerator_as_long(),
                        value.denominator_as_long(),
                    )
                )
            exact_game = tuple(game)
            exact_branch_value = sum(
                (
                    coefficient * value
                    for coefficient, value in zip(
                        branch, exact_game, strict=True
                    )
                ),
                F(0),
            )
            exact_envelope_value = sum(
                (
                    coefficient * value
                    for coefficient, value in zip(
                        envelope, exact_game, strict=True
                    )
                ),
                F(0),
            )
            if exact_game[grand] != 1:
                raise RuntimeError("exact separator normalization failed")
            if any(
                sum(
                    (
                        coefficient * value
                        for coefficient, value in zip(
                            generator, exact_game, strict=True
                        )
                    ),
                    F(0),
                )
                < 0
                for generator in generators
            ):
                raise RuntimeError("exact separator left the game cone")
            if exact_branch_value <= exact_envelope_value:
                raise RuntimeError("exact separator has no positive gap")
            separators.append((exact_game, exact_branch_value))
            print(
                json.dumps(
                    {
                        "exact_separator_z3_fallback": True,
                        "gap": str(
                            exact_branch_value - exact_envelope_value
                        ),
                    }
                ),
                flush=True,
            )
    return separators


def embedded_topology_generators(
    node_ids: list[int],
    n: int,
    protected_pairs: list[dict[str, Any]],
) -> tuple[list[tuple[F, ...]], list[tuple[Any, ...]]]:
    coalition_count = 1 << n
    node_index = {node: index for index, node in enumerate(node_ids)}
    local_generators, local_metadata = exact_generators(n)
    generators: list[tuple[F, ...]] = []
    metadata: list[tuple[Any, ...]] = []

    for node in node_ids:
        offset = node_index[node] * coalition_count
        for generator, row_metadata in zip(
            local_generators, local_metadata, strict=True
        ):
            row = [F(0)] * (len(node_ids) * coalition_count)
            row[offset : offset + coalition_count] = generator
            generators.append(tuple(row))
            metadata.append(("node_game", node, *row_metadata))

    for pair in protected_pairs:
        source = int(pair["source"])
        sink = int(pair["sink"])
        player_mask = int(pair["player_mask"])
        source_offset = node_index[source] * coalition_count
        sink_offset = node_index[sink] * coalition_count
        for coalition in range(coalition_count):
            row = [F(0)] * (len(node_ids) * coalition_count)
            row[source_offset + coalition] = F(-1)
            row[sink_offset + coalition] = F(1)
            if coalition & player_mask == player_mask:
                generators.append(tuple(row))
                metadata.append(
                    (
                        "directed_monotonicity",
                        source,
                        sink,
                        coalition,
                    )
                )
            else:
                generators.append(tuple(row))
                metadata.append(
                    (
                        "protected_invariance",
                        source,
                        sink,
                        coalition,
                        1,
                    )
                )
                generators.append(tuple(-value for value in row))
                metadata.append(
                    (
                        "protected_invariance",
                        source,
                        sink,
                        coalition,
                        -1,
                    )
                )
    return generators, metadata


def extract_envelopes(
    report: dict[str, Any],
    support: dict[str, Any],
    max_denominator: int,
    zero_threshold: float,
) -> dict[int, tuple[F, ...]]:
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    coalition_count = 1 << n
    raw_links: dict[int, dict[int, float]] = defaultdict(dict)
    for active in support["equalities"]:
        row = active["row"]
        if row[0] != "perspective_link":
            continue
        raw_links[int(row[1])][int(row[2])] = float(active["weight"])

    envelopes: dict[int, tuple[F, ...]] = {}
    for terminal in terminals:
        node = int(terminal["node"])
        divergence = tuple(F(value) for value in terminal["scaled_divergence"])
        if len(set(divergence)) <= 2:
            envelopes[node] = coefficient_vector(
                choquet_representation(divergence, n), n
            )
            continue
        if node not in raw_links:
            raise RuntimeError(f"support omits perspective links for node {node}")
        envelope = [F(0)] * coalition_count
        for coalition, value in raw_links[node].items():
            if abs(value) <= zero_threshold:
                continue
            envelope[coalition] = F(value).limit_denominator(max_denominator)
        # Every branch representation evaluates to d_i on the additive
        # unanimity ray u_i(S)=1[i in S].  A zero-gap global envelope must do
        # the same.  Floating duals satisfy these identities only up to solver
        # tolerance, so restore them exactly using the singleton coordinates.
        for player, target in enumerate(divergence):
            singleton = 1 << player
            envelope[singleton] = target - sum(
                (
                    envelope[coalition]
                    for coalition in range(coalition_count)
                    if coalition != singleton and coalition >> player & 1
                ),
                F(0),
            )
        envelopes[node] = tuple(envelope)
    return envelopes


def protected_pairs_from_report(report: dict[str, Any]) -> list[dict[str, Any]]:
    players: dict[tuple[int, int], set[int]] = defaultdict(set)
    for path in report["terminal_path_decomposition"]:
        players[int(path["source"]), int(path["sink"])].add(
            int(path["player"])
        )
    return [
        {
            "source": source,
            "sink": sink,
            "players": sorted(block),
            "player_mask": sum(1 << player for player in block),
        }
        for (source, sink), block in sorted(players.items())
    ]


def generate(
    report_path: Path,
    support_path: Path,
    output: Path,
    max_denominator: int,
    zero_threshold: float,
    topology_projection_epsilon: F | None,
    max_preseeded_memberships_per_node: int | None,
) -> None:
    report = json.loads(report_path.read_text())
    raw_support_payload = json.loads(support_path.read_text())
    support_payload = raw_support_payload
    if "perspective_link_support" in raw_support_payload:
        support_payload = raw_support_payload["perspective_link_support"]
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    coalition_count = 1 << n
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    rounded_envelopes = extract_envelopes(
        report,
        support_payload,
        max_denominator,
        zero_threshold,
    )
    protected_pairs = protected_pairs_from_report(report)
    topology_generators, topology_metadata = embedded_topology_generators(
        node_ids, n, protected_pairs
    )
    terminal_divergences = {
        int(terminal["node"]): tuple(
            F(value) for value in terminal["scaled_divergence"]
        )
        for terminal in terminals
    }
    local_generators, local_metadata = exact_generators(n)
    projection_cuts: list[tuple[int, tuple[F, ...], F]] = []
    projection_branch_memberships: list[
        tuple[int, int, tuple[F, ...]]
    ] = []
    seen_preseeded_memberships = set()
    generated_projection_memberships = set()
    projection_preseeded_cuts = raw_support_payload.get(
        "cuts", raw_support_payload.get("active_cuts", [])
    )
    preseeded_cuts = raw_support_payload.get(
        "active_cuts", projection_preseeded_cuts
    )
    seen_projection_cuts = set()
    for active in projection_preseeded_cuts:
        cut = (
            int(active["node"]),
            tuple(F(value) for value in active["game"]),
            F(active["branch_value"]),
        )
        if cut in seen_projection_cuts:
            continue
        seen_projection_cuts.add(cut)
        projection_cuts.append(cut)
    refinement_checkpoint = output.with_suffix(
        output.suffix + ".refinements.json"
    )
    if refinement_checkpoint.exists():
        checkpoint_payload = json.loads(refinement_checkpoint.read_text())
        if (
            checkpoint_payload.get("source_report") == str(report_path)
            and checkpoint_payload.get("source_float_support")
            == str(support_path)
        ):
            for active in checkpoint_payload.get("generated_cuts", []):
                cut = (
                    int(active["node"]),
                    tuple(F(value) for value in active["game"]),
                    F(active["branch_value"]),
                )
                if cut in seen_projection_cuts:
                    continue
                seen_projection_cuts.add(cut)
                projection_cuts.append(cut)
            for active in checkpoint_payload.get(
                "generated_memberships", []
            ):
                node = int(active["node"])
                branch = int(active["branch"])
                membership_pair = (node, branch)
                if membership_pair in seen_preseeded_memberships:
                    continue
                representations = extreme_representations(
                    list(terminal_divergences[node]), n
                )
                if not 0 <= branch < len(representations):
                    raise RuntimeError(
                        "checkpointed branch index is out of range"
                    )
                projection_branch_memberships.append(
                    (
                        node,
                        branch,
                        coefficient_vector(representations[branch], n),
                    )
                )
                seen_preseeded_memberships.add(membership_pair)
                generated_projection_memberships.add(membership_pair)

    def save_refinement_checkpoint() -> None:
        checkpoint_payload = {
            "status": "n5_exact_projection_refinement_checkpoint",
            "source_report": str(report_path),
            "source_float_support": str(support_path),
            "generated_cuts": [
                {
                    "node": node,
                    "game": [str(value) for value in game],
                    "branch_value": str(branch_value),
                }
                for node, game, branch_value in projection_cuts
                if (node, game, branch_value)
                not in {
                    (
                        int(active["node"]),
                        tuple(F(value) for value in active["game"]),
                        F(active["branch_value"]),
                    )
                    for active in projection_preseeded_cuts
                }
            ],
            "generated_memberships": [
                {"node": node, "branch": branch}
                for node, branch in sorted(
                    generated_projection_memberships
                )
            ],
        }
        pending = refinement_checkpoint.with_suffix(
            refinement_checkpoint.suffix + ".pending"
        )
        pending.write_text(json.dumps(checkpoint_payload, indent=2) + "\n")
        pending.replace(refinement_checkpoint)
    if max_preseeded_memberships_per_node is not None:
        preseeded_cuts = sorted(
            preseeded_cuts,
            key=lambda active: (
                int(active["node"]),
                abs(float(active.get("slack", 0.0))),
                int(active["branch"]),
            ),
        )
    preseeded_count_by_node: dict[int, int] = defaultdict(int)
    for active in preseeded_cuts:
        node = int(active["node"])
        branch = int(active["branch"])
        if (node, branch) in seen_preseeded_memberships:
            continue
        if (
            max_preseeded_memberships_per_node is not None
            and preseeded_count_by_node[node]
            >= max_preseeded_memberships_per_node
        ):
            continue
        seen_preseeded_memberships.add((node, branch))
        preseeded_count_by_node[node] += 1
        representations = extreme_representations(
            list(terminal_divergences[node]), n
        )
        if not 0 <= branch < len(representations):
            raise RuntimeError("preseeded branch index is out of range")
        projection_branch_memberships.append(
            (node, branch, coefficient_vector(representations[branch], n))
        )
    for projection_refinement in range(20):
        projected_topology_support = None
        envelopes = rounded_envelopes
        if topology_projection_epsilon is not None:
            last_reconstruction_error = None
            for topology_objective_salt in range(33):
                try:
                    (
                        envelopes,
                        projected_topology_support,
                    ) = exact_topology_projection(
                        rounded_envelopes,
                        terminal_divergences,
                        node_ids,
                        n,
                        topology_generators,
                        topology_projection_epsilon,
                        projection_cuts,
                        projection_branch_memberships,
                        local_generators,
                        topology_objective_salt,
                    )
                    break
                except TopologyReconstructionError as error:
                    last_reconstruction_error = error
                    print(
                        json.dumps(
                            {
                                "topology_projection_alternate_objective": (
                                    topology_objective_salt + 1
                                )
                            }
                        ),
                        flush=True,
                    )
            else:
                raise RuntimeError(
                    "no topology objective admitted exact reconstruction"
                ) from last_reconstruction_error
        # Recheck the memberships that define the current projection first.
        # These are precisely the historically difficult boundary faces, so a
        # failed refinement should be discovered before replaying hundreds of
        # branches already known to be easy.  On the successful final pass the
        # second loop still adds every remaining mixed branch exactly once.
        representations_by_node = {
            node: extreme_representations(
                list(terminal_divergences[node]), n
            )
            for node in node_ids
            if len(set(terminal_divergences[node])) > 2
        }
        branch_jobs = []
        seen_branch_jobs = set()
        for node, branch, _ in projection_branch_memberships:
            if node not in representations_by_node:
                continue
            job_key = (node, branch)
            if job_key in seen_branch_jobs:
                continue
            branch_jobs.append(
                (node, branch, representations_by_node[node][branch])
            )
            seen_branch_jobs.add(job_key)
        for node in node_ids:
            for branch, representation in enumerate(
                representations_by_node.get(node, [])
            ):
                job_key = (node, branch)
                if job_key in seen_branch_jobs:
                    continue
                branch_jobs.append((node, branch, representation))
                seen_branch_jobs.add(job_key)

        local_certificates = []
        failure = None
        for node, branch, representation in branch_jobs:
            branch_vector = coefficient_vector(representation, n)
            target = tuple(
                envelope_value - branch_value
                for envelope_value, branch_value in zip(
                    envelopes[node], branch_vector, strict=True
                )
            )
            try:
                cone_support = find_certificate(
                    target,
                    local_generators,
                    allow_z3=False,
                    alternate_basis_retries=0,
                )
            except RuntimeError as error:
                separators = exact_separating_games(
                    envelopes[node], branch_vector, local_generators, n
                )
                if separators:
                    failure = (
                        node,
                        branch,
                        branch_vector,
                        error,
                        separators,
                    )
                    break
                try:
                    cone_support = find_certificate(
                        target, local_generators
                    )
                except RuntimeError as recovered_error:
                    failure = (
                        node,
                        branch,
                        branch_vector,
                        recovered_error,
                        None,
                    )
                    break
            beta, terms = representation
            local_certificates.append(
                {
                    "node": node,
                    "branch": branch,
                    "beta": str(beta),
                    "terms": [
                        [coalition, str(weight)]
                        for coalition, weight in terms
                    ],
                    "dual_cone_support": [
                        {
                            "generator": list(local_metadata[index]),
                            "weight": str(weight),
                        }
                        for index, weight in cone_support
                    ],
                }
            )
        if failure is None:
            break
        node, branch, branch_vector, error, separators = failure
        if topology_projection_epsilon is None:
            raise RuntimeError(
                f"local envelope failed at node {node}, branch {branch}"
            ) from error
        if separators is None:
            separators = exact_separating_games(
                envelopes[node], branch_vector, local_generators, n
            )
        if not separators:
            raise RuntimeError(
                f"failed to construct an exact separating cut at node "
                f"{node}, branch {branch}"
            ) from error
        new_cuts = [
            (node, game, branch_value)
            for game, branch_value in separators
            if (node, game, branch_value) not in seen_projection_cuts
        ]
        if not new_cuts:
            raise RuntimeError(
                f"all exact separating cuts repeated at node {node}, branch "
                f"{branch}"
            ) from error
        projection_cuts.extend(new_cuts)
        seen_projection_cuts.update(new_cuts)
        membership_pair = (node, branch)
        if membership_pair not in seen_preseeded_memberships:
            projection_branch_memberships.append(
                (node, branch, branch_vector)
            )
            seen_preseeded_memberships.add(membership_pair)
            generated_projection_memberships.add(membership_pair)
        save_refinement_checkpoint()
        print(
            json.dumps(
                {
                    "projection_refinement": projection_refinement,
                    "failed_node": node,
                    "failed_branch": branch,
                    "new_exact_projection_cut_count": len(new_cuts),
                    "exact_projection_cut_count": len(projection_cuts),
                    "exact_projection_branch_membership_count": len(
                        projection_branch_memberships
                    ),
                }
            ),
            flush=True,
        )
    else:
        raise RuntimeError("exact topology projection refinement did not converge")

    global_target: list[F] = []
    for node in node_ids:
        global_target.extend(-value for value in envelopes[node])
    topology_support = (
        projected_topology_support
        if projected_topology_support is not None
        else find_certificate(tuple(global_target), topology_generators)
    )
    payload = {
        "status": "n5_terminal_perspective_envelope_exact_certificate",
        "n": n,
        "source_report": str(report_path),
        "source_float_support": str(support_path),
        "rounding_max_denominator": max_denominator,
        "rounding_zero_threshold": zero_threshold,
        "topology_projection_epsilon": (
            None
            if topology_projection_epsilon is None
            else str(topology_projection_epsilon)
        ),
        "max_preseeded_memberships_per_node": (
            max_preseeded_memberships_per_node
        ),
        "topology_projection_refinement_count": (
            len(projection_cuts) + len(projection_branch_memberships)
        ),
        "topology_projection_refinement_cuts": [
            {
                "node": node,
                "game": [str(value) for value in game],
                "branch_value": str(branch_value),
            }
            for node, game, branch_value in projection_cuts
        ],
        "topology_projection_branch_memberships": [
            {"node": node, "branch": branch}
            for node, branch, _ in projection_branch_memberships
        ],
        "node_ids": node_ids,
        "terminal_divergences": {
            str(int(terminal["node"])): terminal["scaled_divergence"]
            for terminal in terminals
        },
        "protected_pairs": protected_pairs,
        "envelopes": {
            str(node): [str(value) for value in envelopes[node]]
            for node in node_ids
        },
        "mixed_branch_count": len(local_certificates),
        "local_envelope_certificates": local_certificates,
        "global_topology_certificate": [
            {
                "generator": list(topology_metadata[index]),
                "weight": str(weight),
            }
            for index, weight in topology_support
        ],
        "exact_maximum": "0",
        "interpretation": (
            "Each envelope dominates every extreme lower-expectation branch "
            "on the monotone exact-game cone. The negative sum of the "
            "envelopes lies in the dual cone of the protected topology, so "
            "the terminal objective is nonpositive for arbitrary grand "
            "worths. The zero game attains equality."
        ),
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")


def reconstruct(
    support: list[dict[str, Any]],
    generator_lookup: dict[tuple[Any, ...], int],
    generators: list[tuple[F, ...]],
) -> tuple[F, ...]:
    result = [F(0)] * len(generators[0])
    for active in support:
        metadata = tuple(active["generator"])
        if metadata not in generator_lookup:
            raise RuntimeError(f"unknown generator {metadata}")
        weight = F(active["weight"])
        if weight < 0:
            raise RuntimeError("negative dual-cone coefficient")
        generator = generators[generator_lookup[metadata]]
        for coordinate, coefficient in enumerate(generator):
            result[coordinate] += weight * coefficient
    return tuple(result)


def verify(archive: Path) -> None:
    payload = json.loads(archive.read_text())
    if payload["status"] != "n5_terminal_perspective_envelope_exact_certificate":
        raise RuntimeError("unexpected archive status")
    n = int(payload["n"])
    node_ids = [int(node) for node in payload["node_ids"]]
    envelopes = {
        int(node): tuple(F(value) for value in values)
        for node, values in payload["envelopes"].items()
    }
    local_generators, local_metadata = exact_generators(n)
    local_lookup = {
        tuple(metadata): index
        for index, metadata in enumerate(local_metadata)
    }
    certificate_lookup = {
        (int(entry["node"]), int(entry["branch"])): entry
        for entry in payload["local_envelope_certificates"]
    }
    expected_branch_count = 0
    for node in node_ids:
        divergence = tuple(
            F(value)
            for value in payload["terminal_divergences"][str(node)]
        )
        if len(set(divergence)) <= 2:
            expected = coefficient_vector(
                choquet_representation(divergence, n), n
            )
            if envelopes[node] != expected:
                raise RuntimeError("two-level envelope is not the exact identity")
            continue
        representations = extreme_representations(list(divergence), n)
        expected_branch_count += len(representations)
        for branch, representation in enumerate(representations):
            key = (node, branch)
            if key not in certificate_lookup:
                raise RuntimeError(f"missing local certificate {key}")
            entry = certificate_lookup[key]
            archived_representation = (
                F(entry["beta"]),
                tuple(
                    (int(coalition), F(weight))
                    for coalition, weight in entry["terms"]
                ),
            )
            if archived_representation != representation:
                raise RuntimeError("branch representation mismatch")
            target = tuple(
                envelope_value - branch_value
                for envelope_value, branch_value in zip(
                    envelopes[node],
                    coefficient_vector(representation, n),
                    strict=True,
                )
            )
            if reconstruct(
                entry["dual_cone_support"],
                local_lookup,
                local_generators,
            ) != target:
                raise RuntimeError("local envelope identity failed")
    if len(certificate_lookup) != expected_branch_count:
        raise RuntimeError("duplicate or extra local branch certificate")
    if int(payload["mixed_branch_count"]) != expected_branch_count:
        raise RuntimeError("mixed branch count mismatch")

    topology_generators, topology_metadata = embedded_topology_generators(
        node_ids, n, payload["protected_pairs"]
    )
    topology_lookup = {
        tuple(metadata): index
        for index, metadata in enumerate(topology_metadata)
    }
    global_target: list[F] = []
    for node in node_ids:
        global_target.extend(-value for value in envelopes[node])
    if reconstruct(
        payload["global_topology_certificate"],
        topology_lookup,
        topology_generators,
    ) != tuple(global_target):
        raise RuntimeError("global topology identity failed")
    if payload["exact_maximum"] != "0":
        raise RuntimeError("unexpected exact maximum")
    print(
        "PASS: exact terminal-envelope certificate, "
        f"{expected_branch_count} mixed branches"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--generate-from-report", type=Path)
    parser.add_argument("--link-support", type=Path)
    parser.add_argument("--max-denominator", type=int, default=8)
    parser.add_argument("--zero-threshold", type=float, default=0.02)
    parser.add_argument("--topology-projection-epsilon", type=F)
    parser.add_argument(
        "--max-preseeded-memberships-per-node", type=int, default=2
    )
    args = parser.parse_args()
    if (
        args.max_preseeded_memberships_per_node is not None
        and args.max_preseeded_memberships_per_node < 0
    ):
        parser.error("--max-preseeded-memberships-per-node must be nonnegative")
    if args.generate_from_report is not None:
        if args.link_support is None:
            raise RuntimeError("--link-support is required when generating")
        generate(
            args.generate_from_report,
            args.link_support,
            args.archive,
            args.max_denominator,
            args.zero_threshold,
            args.topology_projection_epsilon,
            args.max_preseeded_memberships_per_node,
        )
    verify(args.archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
