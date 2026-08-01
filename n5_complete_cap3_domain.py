#!/usr/bin/env python3
"""Solve the complete monotone exact five-player integer grid through cap 3."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_boolean_exact_domain import GRAND, N
from n5_facet_search import load_facets
from n5_quaternary_reachable_domain import exact_grand_two_games, induced_edges
from n5_symmetric_family_lp import integer_grid_certificate, symmetric_family_result


class Collector(cp_model.CpSolverSolutionCallback):
    def __init__(self, variables) -> None:
        super().__init__()
        self.variables = variables
        self.games: list[bytes] = []

    def on_solution_callback(self) -> None:
        self.games.append(bytes(self.value(variable) for variable in self.variables))


def enumerate_grand_three(facets) -> list[bytes]:
    model = cp_model.CpModel()
    values = [model.new_int_var(0, 3, f"v_{coalition}") for coalition in range(GRAND + 1)]
    model.add(values[0] == 0)
    model.add(values[GRAND] == 3)
    for coalition in range(GRAND + 1):
        for player in range(N):
            if not coalition >> player & 1:
                model.add(values[coalition] <= values[coalition | (1 << player)])
    for facet in facets:
        model.add(
            sum(
                int(facet[coalition]) * values[coalition]
                for coalition in range(1, GRAND + 1)
                if facet[coalition]
            )
            >= 0
        )
    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_search_workers = 1
    collector = Collector(values)
    status = solver.solve(model, collector)
    if status != cp_model.OPTIMAL:
        raise RuntimeError(f"complete enumeration failed: {solver.status_name(status)}")
    return collector.games


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facets, _types = load_facets()
    facet_matrix = np.asarray(facets, dtype=np.int16)[:, 1:]
    lower = exact_grand_two_games(facet_matrix)
    grand_three = enumerate_grand_three(facets)
    games = sorted(lower) + sorted(grand_three)
    edges = induced_edges(games, 3, True)
    result = symmetric_family_result(games, edges, "highs-ds")
    if not result.success:
        raise RuntimeError(result.message)
    margin = float(result.x[-1])
    exact_verification = integer_grid_certificate(
        games, sorted(edges), result.quotient_variables
    )
    if not exact_verification["verified"]:
        raise RuntimeError(f"integer certificate failed: {exact_verification}")
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if margin < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "scope": "complete monotone exact five-player integer grid with worths in {0,1,2,3}",
        "lower_section_game_count": len(lower),
        "grand_three_game_count": len(grand_three),
        "game_count": len(games),
        "legal_edge_count": len(edges),
        "symmetry_quotient": result.quotient_counts,
        "margin": margin,
        "exact_integer_primal_certificate": exact_verification,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if margin < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
