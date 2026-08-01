#!/usr/bin/env python3
"""Globally optimize fixed flows on legal exact bowtie complexes."""

from __future__ import annotations

import argparse
import json
import random
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


PROJECT_SRC = Path(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)
sys.path.insert(0, str(PROJECT_SRC))

from n5_facet_search import load_facets

from n5_mixed_terminal_milp import lower_dual_vertices


N = 5
GRAND = 31


@dataclass(frozen=True)
class FlowTopology:
    node_count: int
    edges: tuple[tuple[int, int, int], ...]
    flows: tuple[tuple[int, ...], ...]
    divergence: tuple[tuple[int, ...], ...]


def bowtie(
    directions: tuple[tuple[int, int], ...], rng: random.Random
) -> FlowTopology:
    edges = []
    for arm, (first, second) in enumerate(directions):
        left = 1 + 3 * arm
        right = left + 1
        top = left + 2
        edges.extend(
            ((0, left, first), (0, right, second), (left, top, second), (right, top, first))
        )
    flows = []
    divergence = [[0] * N for _ in range(1 + 3 * len(directions))]
    for lower, upper, coalition in edges:
        flow = tuple(
            rng.randint(1, 3) if coalition >> player & 1 else 0
            for player in range(N)
        )
        flows.append(flow)
        for player, weight in enumerate(flow):
            divergence[lower][player] += weight
            divergence[upper][player] -= weight
    return FlowTopology(
        len(divergence), tuple(edges), tuple(flows), tuple(tuple(row) for row in divergence)
    )


def capped_bowtie(
    directions: tuple[tuple[int, int], ...], rng: random.Random
) -> FlowTopology:
    edges = [(0, 1, GRAND)]
    for arm, (first, second) in enumerate(directions):
        left = 2 + 3 * arm
        right = left + 1
        top = left + 2
        edges.extend(
            (
                (1, left, first),
                (1, right, second),
                (left, top, second),
                (right, top, first),
            )
        )
    flows = []
    divergence = [[0] * N for _ in range(2 + 3 * len(directions))]
    for lower, upper, coalition in edges:
        flow = tuple(
            rng.randint(1, 3) if coalition >> player & 1 else 0
            for player in range(N)
        )
        flows.append(flow)
        for player, weight in enumerate(flow):
            divergence[lower][player] += weight
            divergence[upper][player] -= weight
    return FlowTopology(
        len(divergence), tuple(edges), tuple(flows), tuple(tuple(row) for row in divergence)
    )


class LegalExactCone:
    def __init__(
        self,
        topology: FlowTopology,
        minimum_bump: float,
        normalization_node: int = 0,
        fixed_bump: bool = False,
        tangent_player: int | None = None,
    ) -> None:
        self.topology = topology
        self.tangent_player = tangent_player
        self.coalition_count = GRAND + 1
        self.game_variables = topology.node_count * self.coalition_count
        self.rows: list[dict[int, float]] = []
        self.lower: list[float] = []
        self.upper: list[float] = []
        facets, _ = load_facets()

        def add(row, lower=-np.inf, upper=np.inf):
            self.rows.append({column: value for column, value in row.items() if value})
            self.lower.append(lower)
            self.upper.append(upper)

        for node in range(topology.node_count):
            add({self.column(node, 0): 1.0}, 0.0, 0.0)
            if tangent_player is not None:
                add({self.column(node, GRAND): 1.0}, 0.0, 0.0)
            elif node == normalization_node:
                add({self.column(node, GRAND): 1.0}, 1.0, 1.0)
            for facet in facets:
                add(
                    {
                        self.column(node, coalition): -float(coefficient)
                        for coalition, coefficient in enumerate(facet)
                        if coefficient
                    },
                    upper=0.0,
                )
            for coalition in range(GRAND + 1):
                for player in range(N):
                    if coalition >> player & 1:
                        continue
                    if tangent_player is not None and (
                        tangent_player < 0 or player == tangent_player
                    ):
                        continue
                    successor = coalition | (1 << player)
                    add(
                        {
                            self.column(node, coalition): 1.0,
                            self.column(node, successor): -1.0,
                        },
                        upper=0.0,
                    )
        for lower_node, upper_node, changed in topology.edges:
            for coalition in range(GRAND + 1):
                row = {
                    self.column(lower_node, coalition): 1.0,
                    self.column(upper_node, coalition): -1.0,
                }
                if coalition == changed:
                    if fixed_bump:
                        add(row, -minimum_bump, -minimum_bump)
                    else:
                        add(row, upper=-minimum_bump)
                else:
                    add(row, 0.0, 0.0)

    def column(self, node: int, coalition: int) -> int:
        return node * self.coalition_count + coalition


def solve(
    topology: FlowTopology,
    time_limit: float,
    minimum_bump: float,
    tangent_player: int | None = None,
) -> dict[str, object]:
    cone = LegalExactCone(
        topology, minimum_bump, tangent_player=tangent_player
    )
    candidates = [lower_dual_vertices(vector, N) for vector in topology.divergence]
    z_start = cone.game_variables
    binary_starts = []
    variable_count = cone.game_variables + topology.node_count
    for node_candidates in candidates:
        binary_starts.append(variable_count)
        variable_count += len(node_candidates)

    rows = list(cone.rows)
    lower = list(cone.lower)
    upper = list(cone.upper)

    def add(row, low=-np.inf, high=np.inf):
        rows.append({column: value for column, value in row.items() if value})
        lower.append(low)
        upper.append(high)

    for node, (divergence, node_candidates) in enumerate(
        zip(topology.divergence, candidates, strict=True)
    ):
        z = z_start + node
        minimum = float(min(divergence))
        maximum = float(max(divergence))
        for candidate_index, candidate in enumerate(node_candidates):
            constant = float(candidate["constant"])
            coefficients = candidate["coefficients"]
            row = {z: -1.0}
            for coalition, coefficient in coefficients.items():
                row[cone.column(node, coalition)] = float(coefficient)
            add(row, high=-constant)

            if tangent_player is None:
                big_m = max(1.0, maximum - constant)
            else:
                expression_lower = constant - sum(
                    abs(float(coefficient))
                    for coefficient in coefficients.values()
                )
                big_m = maximum - expression_lower
            row = {
                z: 1.0,
                binary_starts[node] + candidate_index: big_m,
            }
            for coalition, coefficient in coefficients.items():
                row[cone.column(node, coalition)] = -float(coefficient)
            add(row, high=big_m + constant)
        add(
            {
                binary_starts[node] + candidate_index: 1.0
                for candidate_index in range(len(node_candidates))
            },
            1.0,
            1.0,
        )

    rr, cc, vv = [], [], []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            rr.append(row_index)
            cc.append(column)
            vv.append(value)
    matrix = coo_matrix((vv, (rr, cc)), shape=(len(rows), variable_count)).tocsr()
    objective = np.zeros(variable_count)
    objective[z_start : z_start + topology.node_count] = -1.0
    variable_lower = np.full(variable_count, -np.inf)
    variable_upper = np.full(variable_count, np.inf)
    variable_lower[: cone.game_variables] = (
        -1.0 if tangent_player is not None else 0.0
    )
    variable_upper[: cone.game_variables] = 1.0
    integrality = np.zeros(variable_count)
    for node, (vector, node_candidates) in enumerate(zip(topology.divergence, candidates)):
        if tangent_player is None:
            variable_lower[z_start + node] = min(vector)
            variable_upper[z_start + node] = max(vector)
        else:
            bound = sum(abs(value) for value in vector)
            variable_lower[z_start + node] = -bound
            variable_upper[z_start + node] = bound
        start = binary_starts[node]
        stop = start + len(node_candidates)
        variable_lower[start:stop] = 0.0
        variable_upper[start:stop] = 1.0
        integrality[start:stop] = 1
    result = milp(
        objective,
        integrality=integrality,
        bounds=Bounds(variable_lower, variable_upper),
        constraints=LinearConstraint(matrix, np.asarray(lower), np.asarray(upper)),
        options={"time_limit": time_limit, "mip_rel_gap": 0.0},
    )
    payload: dict[str, object] = {
        "status": "optimal" if result.success else "incomplete",
        "message": result.message,
        "objective": None if result.fun is None else -float(result.fun),
        "mip_gap": getattr(result, "mip_gap", None),
        "candidate_counts": [len(value) for value in candidates],
        "edges": [list(edge) for edge in topology.edges],
        "flows": [list(flow) for flow in topology.flows],
        "divergence": [list(vector) for vector in topology.divergence],
        "tangent_player": tangent_player,
    }
    if result.x is not None:
        payload["node_lower_expectations"] = [
            float(result.x[z_start + node]) for node in range(topology.node_count)
        ]
        payload["games"] = [
            [
                float(result.x[cone.column(node, coalition)])
                for coalition in range(GRAND + 1)
            ]
            for node in range(topology.node_count)
        ]
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=10)
    parser.add_argument("--arms", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--minimum-bump", type=float, default=0.001)
    parser.add_argument("--tangent-player", type=int, choices=range(-1, N))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pool = [coalition for coalition in range(1, GRAND) if 2 <= coalition.bit_count() <= 4]
    rows = []
    best = None
    for case in range(args.cases):
        directions = tuple(tuple(rng.sample(pool, 2)) for _ in range(args.arms))
        topology = bowtie(directions, rng)
        payload = solve(
            topology,
            args.time_limit,
            args.minimum_bump,
            args.tangent_player,
        )
        row = {"case": case, "directions": directions, **payload}
        rows.append(row)
        if payload["objective"] is not None and (
            best is None or payload["objective"] > best["objective"]
        ):
            best = row
        print(
            json.dumps(
                {
                    "case": case,
                    "directions": directions,
                    "status": payload["status"],
                    "objective": payload["objective"],
                    "candidate_counts": payload["candidate_counts"],
                }
            ),
            flush=True,
        )
        if payload["objective"] is not None and payload["objective"] > 1e-8:
            break
    output = {
        "status": (
            "positive_flow_found"
            if best is not None and best["objective"] > 1e-8
            else "none_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "best": best,
        "rows": rows,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
