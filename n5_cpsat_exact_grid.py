#!/usr/bin/env python3
"""Enumerate bounded five-player exact games with CP-SAT."""

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


class Counter(cp_model.CpSolverSolutionCallback):
    def __init__(self, variables, progress: int) -> None:
        super().__init__()
        self.variables = variables
        self.progress = progress
        self.count = 0
        self.examples = []

    def on_solution_callback(self) -> None:
        self.count += 1
        if len(self.examples) < 10:
            self.examples.append([self.value(variable) for variable in self.variables])
        if self.count % self.progress == 0:
            print(
                json.dumps(
                    {
                        "solutions": self.count,
                        "wall_time": self.wall_time,
                        "branches": self.num_branches,
                    }
                ),
                flush=True,
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cap", type=int, default=3)
    parser.add_argument("--seconds", type=float, default=60.0)
    parser.add_argument("--progress", type=int, default=100_000)
    parser.add_argument("--order-singletons", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    facets, _types = load_facets()
    model = cp_model.CpModel()
    values = [model.new_int_var(0, args.cap, f"v_{coalition}") for coalition in range(GRAND + 1)]
    model.add(values[0] == 0)
    model.add(values[GRAND] == args.cap)
    for coalition in range(GRAND + 1):
        for player in range(N):
            if coalition >> player & 1:
                continue
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
    if args.order_singletons:
        for player in range(N - 1):
            model.add(values[1 << player] <= values[1 << (player + 1)])
    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.max_time_in_seconds = args.seconds
    solver.parameters.num_search_workers = 1
    callback = Counter(values, args.progress)
    status = solver.solve(model, callback)
    payload = {
        "status": solver.status_name(status),
        "cap": args.cap,
        "ordered_singletons": args.order_singletons,
        "solution_count": callback.count,
        "wall_time": solver.wall_time,
        "branches": solver.num_branches,
        "conflicts": solver.num_conflicts,
        "examples": callback.examples,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
