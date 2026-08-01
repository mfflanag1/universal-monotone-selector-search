#!/usr/bin/env python3
"""Adversarially reoptimize exact games on a fixed equivariant bump topology."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_equivariant_boundary_closure import load_quotient
from n5_equivariant_quotient_lp import solve_quotient
from n5_facet_search import load_facets
from n5_intrinsic_permutation_union import permute_game, permute_mask


def sparse(rows, columns):
    row_indices = []
    column_indices = []
    values = []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                row_indices.append(row_index)
                column_indices.append(column)
                values.append(float(value))
    return coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(rows), columns),
    ).tocsr()


def inverse_permutation(permutation):
    inverse = [0] * len(permutation)
    for source, target in enumerate(permutation):
        inverse[target] = source
    return tuple(inverse)


def inverse_mask(mask, permutation):
    return permute_mask(mask, inverse_permutation(permutation))


def raw_witnesses(games, needed_edges):
    permutations = list(itertools.permutations(range(5)))
    orbit_games = []
    orbit_maps = []
    orbit_index = {}
    for game_index, game in enumerate(games):
        for permutation in permutations:
            image = permute_game(game, permutation)
            if image in orbit_index:
                continue
            orbit_index[image] = len(orbit_games)
            orbit_games.append(image)
            orbit_maps.append((game_index, permutation))
    raw_edges = set()
    add_induced_edges_fast(orbit_games, raw_edges)
    witnesses = {}
    for lower, upper, coalition in raw_edges:
        lower_game, lower_permutation = orbit_maps[lower]
        upper_game, upper_permutation = orbit_maps[upper]
        lower_inverse = inverse_permutation(lower_permutation)
        upper_inverse = inverse_permutation(upper_permutation)
        represented = []
        for player in range(5):
            if not coalition >> player & 1:
                continue
            edge = (
                lower_game,
                lower_inverse[player],
                upper_game,
                upper_inverse[player],
            )
            if edge in needed_edges:
                represented.append(edge)
        if not represented:
            continue
        witness = (
            lower_game,
            lower_permutation,
            upper_game,
            upper_permutation,
            coalition,
            orbit_games[upper][coalition] - orbit_games[lower][coalition],
        )
        for edge in represented:
            witnesses.setdefault(edge, witness)
    missing = needed_edges.difference(witnesses)
    if missing:
        raise RuntimeError(f"missing {len(missing)} raw quotient witnesses")
    return list(dict.fromkeys(witnesses.values())), len(orbit_games), len(raw_edges)


def game_model(
    games,
    stabilizers,
    witnesses,
    variable_bumps,
    min_bump,
    tangent_all=False,
):
    game_count = len(games)
    columns = game_count * 32

    def column(game, coalition):
        return game * 32 + coalition

    facets = load_facets()[0]
    ub_rows = []
    ub_rhs = []
    eq_rows = []
    eq_rhs = []
    for game in range(game_count):
        eq_rows.append({column(game, 0): 1.0})
        eq_rhs.append(0.0)
        eq_rows.append({column(game, 31): 1.0})
        eq_rhs.append(0.0 if tangent_all else 1.0)
        for facet in facets:
            ub_rows.append(
                {
                    column(game, coalition): -coefficient
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                }
            )
            ub_rhs.append(0.0)
        if not tangent_all:
            for coalition in range(32):
                for player in range(5):
                    if coalition >> player & 1:
                        continue
                    ub_rows.append(
                        {
                            column(game, coalition): 1.0,
                            column(game, coalition | (1 << player)): -1.0,
                        }
                    )
                    ub_rhs.append(0.0)
        for permutation in stabilizers[game]:
            for coalition in range(32):
                image = permute_mask(coalition, permutation)
                if image <= coalition:
                    continue
                eq_rows.append(
                    {
                        column(game, coalition): 1.0,
                        column(game, image): -1.0,
                    }
                )
                eq_rhs.append(0.0)
    for (
        lower,
        lower_permutation,
        upper,
        upper_permutation,
        changed,
        delta,
    ) in witnesses:
        for coalition in range(32):
            lower_mask = inverse_mask(coalition, lower_permutation)
            upper_mask = inverse_mask(coalition, upper_permutation)
            row = {
                column(upper, upper_mask): 1.0,
                column(lower, lower_mask): -1.0,
            }
            if variable_bumps and coalition == changed:
                ub_rows.append({column: -value for column, value in row.items()})
                ub_rhs.append(-min_bump)
            else:
                eq_rows.append(row)
                eq_rhs.append(float(delta) if coalition == changed else 0.0)
    return (
        sparse(ub_rows, columns),
        np.asarray(ub_rhs),
        sparse(eq_rows, columns),
        np.asarray(eq_rhs),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--variable-bumps", action="store_true")
    parser.add_argument("--min-bump", type=float, default=1e-8)
    parser.add_argument(
        "--selector-methods", default="highs-ipm,highs-ds"
    )
    parser.add_argument("--objective-scale", type=float, default=1e8)
    parser.add_argument("--tangent-all", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    games, edges, stabilizers = load_quotient(args.input)
    witnesses, orbit_count, raw_edge_count = raw_witnesses(games, edges)
    print(
        json.dumps(
            {
                "phase": "witnesses",
                "games": len(games),
                "quotient_edges": len(edges),
                "witnesses": len(witnesses),
                "orbit_games": orbit_count,
                "raw_edges": raw_edge_count,
            }
        ),
        flush=True,
    )
    a_ub, b_ub, a_eq, b_eq = game_model(
        games,
        stabilizers,
        witnesses,
        args.variable_bumps,
        args.min_bump,
        args.tangent_all,
    )
    values = np.asarray([[float(value) for value in game] for game in games])
    if args.tangent_all:
        values = values - np.rint(values)
    best_margin = float("inf")
    best_values = values.copy()
    trace = []
    selector_methods = [
        method.strip()
        for method in args.selector_methods.split(",")
        if method.strip()
    ]
    for iteration in range(args.iterations):
        selector_results = []
        for method in selector_methods:
            margin, selector_result = solve_quotient(
                [tuple(row) for row in values],
                edges,
                stabilizers,
                5,
                method=method,
                return_result=True,
            )
            selector_results.append((margin, method, selector_result))
        margin = min(result[0] for result in selector_results)
        trace.append(margin)
        if margin < best_margin:
            best_margin = margin
            best_values = values.copy()
        print(
            json.dumps(
                {
                    "iteration": iteration,
                    "margin": margin,
                    "best_margin": best_margin,
                }
            ),
            flush=True,
        )
        candidates = []
        for _, method, selector_result in selector_results:
            objective = np.zeros(len(games) * 32)
            marginals = np.asarray(selector_result.ineqlin.marginals)
            for game in range(len(games)):
                for coalition in range(1, 31):
                    objective[game * 32 + coalition] = marginals[
                        game * 30 + coalition - 1
                    ]
            game_result = linprog(
                args.objective_scale * objective,
                A_ub=a_ub,
                b_ub=b_ub,
                A_eq=a_eq,
                b_eq=b_eq,
                bounds=(
                    [(-1.0, 1.0)] * (len(games) * 32)
                    if args.tangent_all
                    else [(None, None)] * (len(games) * 32)
                ),
                method="highs-ds",
                options={
                    "primal_feasibility_tolerance": 1e-10,
                    "dual_feasibility_tolerance": 1e-10,
                },
            )
            if not game_result.success:
                raise RuntimeError(game_result.message)
            candidate = game_result.x.reshape(len(games), 32)
            residual = max(
                float(np.max(a_ub @ game_result.x - b_ub)),
                float(np.max(np.abs(a_eq @ game_result.x - b_eq))),
            )
            if residual > 1e-8:
                raise RuntimeError(f"game LP residual {residual}")
            candidate_margins = [
                solve_quotient(
                    [tuple(row) for row in candidate],
                    edges,
                    stabilizers,
                    5,
                    method=selector_method,
                )
                for selector_method in selector_methods
            ]
            candidates.append(
                (min(candidate_margins), method, candidate, residual)
            )
        updated_margin, driving_method, updated, residual = min(
            candidates, key=lambda candidate: candidate[0]
        )
        print(
            json.dumps(
                {
                    "iteration": iteration,
                    "candidate_margin": updated_margin,
                    "driving_method": driving_method,
                    "game_residual": residual,
                }
            ),
            flush=True,
        )
        if updated_margin >= margin - 1e-12:
            trace.append(updated_margin)
            if updated_margin < best_margin:
                best_margin = updated_margin
                best_values = updated.copy()
            break
        values = updated
    output = {
        "status": (
            "equivariant_obstruction_candidate"
            if best_margin < -1e-8
            else "no_equivariant_obstruction_found"
        ),
        "source": str(args.input),
        "game_count": len(games),
        "quotient_edge_count": len(edges),
        "witness_count": len(witnesses),
        "variable_bumps": args.variable_bumps,
        "min_bump": args.min_bump,
        "selector_methods": selector_methods,
        "objective_scale": args.objective_scale,
        "tangent_all": args.tangent_all,
        "game_inequalities": a_ub.shape[0],
        "game_equalities": a_eq.shape[0],
        "trace": trace,
        "best_margin": best_margin,
        "best_games_float": best_values.tolist(),
        "quotient_edges": [list(edge) for edge in sorted(edges)],
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: value for key, value in output.items() if key != "best_games_float"}, indent=2))
    return 1 if best_margin < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
