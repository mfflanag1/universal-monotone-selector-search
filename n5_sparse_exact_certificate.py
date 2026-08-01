#!/usr/bin/env python3
"""Create an exact n=5 family archive from sparse primal-dual certificates."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, csr_matrix

from family_search import Family, common_core_gap_exact, protected_path_width_bound
from n5_facet_search import load_facets, topology_depths


F = Fraction
N = 5
GRAND = (1 << N) - 1


@dataclass(frozen=True)
class SparseRow:
    entries: tuple[tuple[int, F], ...]
    rhs: F
    metadata: tuple[Any, ...]


def members(coalition: int) -> tuple[int, ...]:
    return tuple(player for player in range(N) if coalition >> player & 1)


def load_record(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    record = payload.get("record")
    if record is not None:
        return record
    search = payload.get("search")
    if search is not None and (
        "best_games_exact" in search or "best_games_float" in search
    ):
        return search
    family = payload.get("family")
    if family is not None:
        return {
            "best_games_exact": family["games"],
            "edges": [
                [
                    int(edge["lower"]),
                    int(edge["upper"]),
                    int(edge["coalition"]),
                ]
                for edge in family["edges"]
            ],
            "source_status": payload.get("status"),
            "source": str(path),
        }
    return payload


def rational_games(record: dict[str, Any], max_denominator: int) -> list[tuple[F, ...]]:
    raw = record.get("best_games_exact", record.get("best_games_float"))
    if raw is None:
        raise ValueError("input has no game family")
    return [
        tuple(F(value).limit_denominator(max_denominator) for value in game)
        for game in raw
    ]


def sparse_constraints(
    games: Sequence[Sequence[F]],
    edges: Sequence[tuple[int, int, int]],
    box: bool,
) -> tuple[list[SparseRow], list[SparseRow]]:
    inequalities: list[SparseRow] = []
    equalities: list[SparseRow] = []
    for game_index, game in enumerate(games):
        equalities.append(
            SparseRow(
                tuple((game_index * N + player, F(1)) for player in range(N)),
                F(game[GRAND]),
                ("efficiency", game_index),
            )
        )
        if box:
            for player in range(N):
                variable = game_index * N + player
                inequalities.append(
                    SparseRow(
                        ((variable, F(-1)),),
                        -F(game[1 << player]),
                        ("box_lower", game_index, player),
                    )
                )
                inequalities.append(
                    SparseRow(
                        ((variable, F(1)),),
                        F(game[GRAND] - game[GRAND ^ (1 << player)]),
                        ("box_upper", game_index, player),
                    )
                )
        else:
            for coalition in range(1, GRAND):
                inequalities.append(
                    SparseRow(
                        tuple(
                            (game_index * N + player, F(-1))
                            for player in members(coalition)
                        ),
                        -F(game[coalition]),
                        ("core", game_index, coalition),
                    )
                )
    margin_index = len(games) * N
    for lower, upper, coalition in edges:
        for player in members(coalition):
            inequalities.append(
                SparseRow(
                    (
                        (lower * N + player, F(1)),
                        (upper * N + player, F(-1)),
                        (margin_index, F(1)),
                    ),
                    F(0),
                    ("monotonicity", lower, upper, coalition, player),
                )
            )
    return inequalities, equalities


def float_matrix(rows: Sequence[SparseRow], variable_count: int) -> csr_matrix:
    row_indices: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for row_index, row in enumerate(rows):
        for column, value in row.entries:
            row_indices.append(row_index)
            columns.append(column)
            values.append(float(value))
    return coo_matrix(
        (values, (row_indices, columns)),
        shape=(len(rows), variable_count),
    ).tocsr()


def exact_lhs(row: SparseRow, point: Sequence[F]) -> F:
    return sum((value * point[column] for column, value in row.entries), F(0))


def reconstruct_primal(
    solution: Sequence[float],
    games: Sequence[Sequence[F]],
    inequalities: Sequence[SparseRow],
    equalities: Sequence[SparseRow],
    denominators: Sequence[int],
    fixed_margin: F | None = None,
) -> tuple[F, list[tuple[F, ...]], list[F]] | None:
    margin_index = len(games) * N
    for max_denominator in denominators:
        point = [F(value).limit_denominator(max_denominator) for value in solution]
        if fixed_margin is not None:
            point[margin_index] = fixed_margin
        for game_index, game in enumerate(games):
            start = game_index * N
            point[start + N - 1] = F(game[GRAND]) - sum(
                point[start : start + N - 1], F(0)
            )
        if any(exact_lhs(row, point) > row.rhs for row in inequalities):
            continue
        if any(exact_lhs(row, point) != row.rhs for row in equalities):
            continue
        allocations = [
            tuple(point[index * N : (index + 1) * N])
            for index in range(len(games))
        ]
        return point[margin_index], allocations, point
    return None


def repair_primal_on_optimal_face(
    games: Sequence[Sequence[F]],
    inequalities: Sequence[SparseRow],
    equalities: Sequence[SparseRow],
    denominators: Sequence[int],
    a_ub: csr_matrix,
    b_ub: np.ndarray,
    a_eq: csr_matrix,
    b_eq: np.ndarray,
    margin: F,
) -> tuple[F, list[tuple[F, ...]], list[F]] | None:
    variable_count = len(games) * N + 1
    rng = np.random.default_rng(20260731)
    bounds = [(None, None)] * (variable_count - 1) + [
        (float(margin), float(margin))
    ]
    for _ in range(6):
        objective = rng.integers(-7, 8, size=variable_count).astype(float)
        objective[-1] = 0.0
        result = linprog(
            objective,
            A_ub=a_ub,
            b_ub=b_ub,
            A_eq=a_eq,
            b_eq=b_eq,
            bounds=bounds,
            method="highs",
        )
        if not result.success:
            continue
        primal = reconstruct_primal(
            result.x,
            games,
            inequalities,
            equalities,
            denominators,
            fixed_margin=margin,
        )
        if primal is not None:
            return primal
    return None


def reconstruct_dual(
    inequality_marginals: Sequence[float],
    equality_marginals: Sequence[float],
    inequalities: Sequence[SparseRow],
    equalities: Sequence[SparseRow],
    variable_count: int,
    margin: F,
    denominators: Sequence[int],
) -> tuple[list[F], list[F], dict[str, Any]] | None:
    active_inequalities = [
        index for index, value in enumerate(inequality_marginals) if abs(value) > 1e-9
    ]
    active_equalities = [
        index for index, value in enumerate(equality_marginals) if abs(value) > 1e-9
    ]
    objective = [F(0)] * variable_count
    objective[-1] = F(-1)
    for max_denominator in denominators:
        inequality_weights = [F(0)] * len(inequalities)
        equality_weights = [F(0)] * len(equalities)
        for index in active_inequalities:
            inequality_weights[index] = F(
                inequality_marginals[index]
            ).limit_denominator(max_denominator)
        for index in active_equalities:
            equality_weights[index] = F(
                equality_marginals[index]
            ).limit_denominator(max_denominator)
        if any(weight > 0 for weight in inequality_weights):
            continue
        stationarity = [F(0)] * variable_count
        for index in active_inequalities:
            weight = inequality_weights[index]
            for column, value in inequalities[index].entries:
                stationarity[column] += weight * value
        for index in active_equalities:
            weight = equality_weights[index]
            for column, value in equalities[index].entries:
                stationarity[column] += weight * value
        if stationarity != objective:
            continue
        dual_objective = sum(
            (inequality_weights[index] * inequalities[index].rhs for index in active_inequalities),
            F(0),
        ) + sum(
            (equality_weights[index] * equalities[index].rhs for index in active_equalities),
            F(0),
        )
        if dual_objective != -margin:
            continue
        certificate = {
            "objective": str(dual_objective),
            "active_inequalities": [
                {
                    "row": list(inequalities[index].metadata),
                    "weight": str(inequality_weights[index]),
                }
                for index in active_inequalities
            ],
            "active_equalities": [
                {
                    "game": index,
                    "weight": str(equality_weights[index]),
                }
                for index in active_equalities
            ],
        }
        return inequality_weights, equality_weights, certificate
    return None


def solve_and_certify(
    games: Sequence[Sequence[F]],
    edges: Sequence[tuple[int, int, int]],
    box: bool,
    denominators: Sequence[int],
) -> tuple[F, list[tuple[F, ...]], dict[str, Any]]:
    inequalities, equalities = sparse_constraints(games, edges, box)
    variable_count = len(games) * N + 1
    objective = np.zeros(variable_count)
    objective[-1] = -1.0
    a_ub = float_matrix(inequalities, variable_count)
    b_ub = np.asarray([float(row.rhs) for row in inequalities])
    a_eq = float_matrix(equalities, variable_count)
    b_eq = np.asarray([float(row.rhs) for row in equalities])
    result = linprog(
        objective,
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=a_eq,
        b_eq=b_eq,
        bounds=[(None, None)] * variable_count,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    primal = reconstruct_primal(
        result.x, games, inequalities, equalities, denominators
    )
    if primal is None:
        margin = F(result.x[-1]).limit_denominator(max(denominators))
        primal = repair_primal_on_optimal_face(
            games,
            inequalities,
            equalities,
            denominators,
            a_ub,
            b_ub,
            a_eq,
            b_eq,
            margin,
        )
    if primal is None:
        raise RuntimeError(
            "failed to reconstruct or repair an exact sparse primal"
        )
    margin, allocations, _ = primal
    dual = reconstruct_dual(
        result.ineqlin.marginals,
        result.eqlin.marginals,
        inequalities,
        equalities,
        variable_count,
        margin,
        denominators,
    )
    if dual is None:
        raise RuntimeError("failed to reconstruct an exact sparse dual")
    return margin, allocations, dual[2]


def verify_game_family(
    games: Sequence[Sequence[F]], edges: Sequence[tuple[int, int, int]]
) -> None:
    if any(len(game) != GRAND + 1 for game in games):
        raise RuntimeError("a game has the wrong coalition count")
    if len({game[GRAND] for game in games}) != 1:
        raise RuntimeError("games do not share a grand worth")
    for lower, upper, coalition in edges:
        differences = [
            games[upper][target] - games[lower][target]
            for target in range(GRAND + 1)
        ]
        if differences[coalition] <= 0 or any(
            value != 0
            for target, value in enumerate(differences)
            if target != coalition
        ):
            raise RuntimeError("an edge is not an exact positive one-coordinate bump")
    facets, _ = load_facets()
    for game in games:
        if any(
            game[coalition] > game[coalition | (1 << player)]
            for coalition in range(GRAND + 1)
            for player in range(N)
            if not coalition >> player & 1
        ):
            raise RuntimeError("a game is not monotone")
        if any(
            sum((F(coefficient) * game[index] for index, coefficient in enumerate(facet)), F(0)) < 0
            for facet in facets
        ):
            raise RuntimeError("a game violates an exact-cone facet")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shape", required=True)
    parser.add_argument("--game-max-denominator", type=int, default=1_000_000)
    parser.add_argument(
        "--certificate-denominators",
        default="1000,10000,100000,1000000,10000000,100000000",
    )
    args = parser.parse_args()

    record = load_record(args.input)
    games = rational_games(record, args.game_max_denominator)
    edges = [tuple(int(value) for value in edge) for edge in record["edges"]]
    verify_game_family(games, edges)
    denominators = [
        int(value) for value in args.certificate_denominators.split(",")
    ]
    margin, allocations, dual = solve_and_certify(
        games, edges, False, denominators
    )
    box_margin, box_allocations, box_dual = solve_and_certify(
        games, edges, True, denominators
    )
    common_gap = common_core_gap_exact(games, N)
    family = Family(
        games=list(games),
        edges=list(edges),
        depths=topology_depths(len(games), edges),
        shape=args.shape,
    )
    path_bound = protected_path_width_bound(family, N)
    archive = {
        "status": "exact_family_certified_sparse",
        "n": N,
        "node_count": len(games),
        "edge_count": len(edges),
        "rational_reconstruction_max_denominator": max(denominators),
        "complete_exact_cone_facet_count": len(load_facets()[0]),
        "float_source": str(args.input),
        "search": {
            key: value
            for key, value in record.items()
            if key not in ("best_games_exact", "best_games_float", "traces")
        },
        "common_core_budget_gap_exact": str(common_gap),
        "max_min_monotonicity_margin_exact": str(margin),
        "box_max_min_monotonicity_margin_exact": str(box_margin),
        "non_atomic_facet_tax_exact": str(box_margin - margin),
        "protected_path_width_bound": (
            {
                "bound": str(path_bound[0]),
                "player": path_bound[1] + 1,
                "path": path_bound[2],
            }
            if path_bound is not None
            else None
        ),
        "family": {
            "shape": args.shape,
            "games": [[str(value) for value in game] for game in games],
            "edges": [
                {
                    "lower": lower,
                    "upper": upper,
                    "coalition": coalition,
                    "members": [player + 1 for player in members(coalition)],
                    "delta": str(games[upper][coalition] - games[lower][coalition]),
                }
                for lower, upper, coalition in edges
            ],
            "depths": family.depths,
        },
        "exact_allocations": [
            [str(value) for value in allocation] for allocation in allocations
        ],
        "exact_margin_dual_certificate": dual,
        "exact_box_allocations": [
            [str(value) for value in allocation] for allocation in box_allocations
        ],
        "exact_box_margin_dual_certificate": box_dual,
    }
    args.output.write_text(json.dumps(archive, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": archive["status"],
                "node_count": archive["node_count"],
                "edge_count": archive["edge_count"],
                "common_core_budget_gap_exact": archive[
                    "common_core_budget_gap_exact"
                ],
                "max_min_monotonicity_margin_exact": archive[
                    "max_min_monotonicity_margin_exact"
                ],
                "box_max_min_monotonicity_margin_exact": archive[
                    "box_max_min_monotonicity_margin_exact"
                ],
                "non_atomic_facet_tax_exact": archive[
                    "non_atomic_facet_tax_exact"
                ],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
