#!/usr/bin/env python3
"""Adversarial selector-margin search on legal exact bowtie complexes."""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from family_search import family_constraints
from n5_actual_bowtie_milp import GRAND, LegalExactCone, bowtie, capped_bowtie
from n5_actual_bowtie_search import matrices
from n5_actual_bowtie_milp import FlowTopology
from n5_random_exact_hourglass import hourglass
from n5_topology_batch import punctured_cube_edges


N = 5


def random_cube_edges(directions, rng, state_probability):
    states = [
        state
        for state in range(1 << len(directions))
        if state == 0 or rng.random() < state_probability
    ]
    node = {state: index for index, state in enumerate(states)}
    edges = []
    for state in states:
        for direction, coalition in enumerate(directions):
            upper = state | (1 << direction)
            if upper != state and upper in node:
                edges.append((node[state], node[upper], coalition))
    return edges


def random_lattice_edges(directions, levels, rng, state_probability):
    states = [
        state
        for state in itertools.product(range(levels), repeat=len(directions))
        if not any(state) or rng.random() < state_probability
    ]
    node = {state: index for index, state in enumerate(states)}
    edges = []
    for state in states:
        for direction, coalition in enumerate(directions):
            if state[direction] + 1 >= levels:
                continue
            upper = list(state)
            upper[direction] += 1
            upper = tuple(upper)
            if upper in node:
                edges.append((node[state], node[upper], coalition))
    return edges


def connected_lattice_edges(directions, levels, rng, target_states):
    zero = (0,) * len(directions)
    states = {zero}
    frontier = set()

    def add_frontier(state):
        for direction in range(len(directions)):
            for step in (-1, 1):
                value = state[direction] + step
                if value < 0 or value >= levels:
                    continue
                neighbor = list(state)
                neighbor[direction] = value
                neighbor = tuple(neighbor)
                if neighbor not in states:
                    frontier.add(neighbor)

    add_frontier(zero)
    while frontier and len(states) < target_states:
        state = rng.choice(tuple(frontier))
        frontier.remove(state)
        states.add(state)
        add_frontier(state)
    ordered = sorted(states)
    node = {state: index for index, state in enumerate(ordered)}
    edges = []
    for state in ordered:
        for direction, coalition in enumerate(directions):
            if state[direction] + 1 >= levels:
                continue
            upper = list(state)
            upper[direction] += 1
            upper = tuple(upper)
            if upper in node:
                edges.append((node[state], node[upper], coalition))
    return edges


def solve_minimum(cone, matrix_data, objective):
    ub, ub_rhs, eq, eq_rhs = matrix_data
    result = linprog(
        objective,
        A_ub=ub,
        b_ub=ub_rhs,
        A_eq=eq,
        b_eq=eq_rhs,
        bounds=[
            (-1.0, 1.0) if cone.tangent_player is not None else (0.0, 1.0)
        ]
        * cone.game_variables,
        method="highs",
    )
    return result.x if result.success else None


def games_from_point(cone, point):
    return [
        [
            float(point[cone.column(node, coalition)])
            for coalition in range(GRAND + 1)
        ]
        for node in range(cone.topology.node_count)
    ]


def margin_dual(games, edges):
    inequalities, rhs, equalities, eq_rhs, metadata = family_constraints(
        games, edges, N
    )
    variable_count = len(games) * N + 1
    a_ub = []
    for row, tag in zip(inequalities, metadata, strict=True):
        a_ub.append(
            [float(value) for value in row]
            + [1.0 if tag[0] == "monotonicity" else 0.0]
        )
    a_eq = [
        [float(value) for value in row] + [0.0] for row in equalities
    ]
    objective = np.zeros(variable_count)
    objective[-1] = -1.0
    result = linprog(
        objective,
        A_ub=np.asarray(a_ub),
        b_ub=np.asarray(rhs, dtype=float),
        A_eq=np.asarray(a_eq),
        b_eq=np.asarray(eq_rhs, dtype=float),
        bounds=[(None, None)] * variable_count,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return (
        float(result.x[-1]),
        result.ineqlin.marginals,
        result.eqlin.marginals,
        metadata,
    )


def fixed_dual_objective(cone, inequality_weights, equality_weights, metadata):
    objective = np.zeros(cone.game_variables)
    for weight, tag in zip(inequality_weights, metadata, strict=True):
        if weight and tag[0] == "core":
            objective[cone.column(int(tag[1]), int(tag[2]))] += weight
    for node, weight in enumerate(equality_weights):
        objective[cone.column(node, GRAND)] -= weight
    return objective


def optimize(
    topology,
    bump,
    starts,
    iterations,
    rng,
    fixed_bump=False,
    tangent_player=None,
):
    normalization_node = 1 if topology.edges[0][2] == GRAND else 0
    cone = LegalExactCone(
        topology,
        bump,
        normalization_node,
        fixed_bump=fixed_bump,
        tangent_player=tangent_player,
    )
    matrix_data = matrices(cone)
    best = None
    for _start in range(starts):
        random_objective = np.asarray(
            [rng.uniform(-1.0, 1.0) for _ in range(cone.game_variables)]
        )
        point = solve_minimum(cone, matrix_data, random_objective)
        if point is None:
            continue
        for _iteration in range(iterations):
            games = games_from_point(cone, point)
            margin, inequality_weights, equality_weights, metadata = margin_dual(
                games, list(topology.edges)
            )
            if best is None or margin < best["margin"]:
                best = {"margin": margin, "games": games}
            objective = fixed_dual_objective(
                cone, inequality_weights, equality_weights, metadata
            )
            updated = solve_minimum(cone, matrix_data, objective)
            if updated is None:
                break
            updated_games = games_from_point(cone, updated)
            updated_margin = margin_dual(updated_games, list(topology.edges))[0]
            point = updated
            if updated_margin >= margin - 1e-9:
                if updated_margin < best["margin"]:
                    best = {"margin": updated_margin, "games": updated_games}
                break
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=100)
    parser.add_argument("--arms", type=int, default=2)
    parser.add_argument("--starts", type=int, default=5)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--minimum-bump", type=float, default=0.001)
    parser.add_argument("--aggregate-cap", action="store_true")
    parser.add_argument("--hourglass", action="store_true")
    parser.add_argument("--punctured-cube", action="store_true")
    parser.add_argument("--random-cube", action="store_true")
    parser.add_argument("--random-lattice", action="store_true")
    parser.add_argument("--connected-lattice", action="store_true")
    parser.add_argument("--levels", type=int, default=3)
    parser.add_argument("--target-states", type=int, default=40)
    parser.add_argument("--state-probability", type=float, default=0.7)
    parser.add_argument("--fixed-bump", action="store_true")
    parser.add_argument("--tangent-player", type=int, choices=range(-1, N))
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    pool = [mask for mask in range(1, GRAND) if 2 <= mask.bit_count() <= 4]
    rows = []
    best = None
    for case in range(args.cases):
        if args.connected_lattice:
            directions = tuple(rng.sample(pool, args.arms))
            topology_edges = tuple(
                connected_lattice_edges(
                    directions,
                    args.levels,
                    rng,
                    args.target_states,
                )
            )
            topology = FlowTopology(
                1 + max(max(edge[0], edge[1]) for edge in topology_edges),
                topology_edges,
                tuple(),
                tuple(),
            )
        elif args.random_lattice:
            directions = tuple(rng.sample(pool, args.arms))
            topology_edges = tuple(
                random_lattice_edges(
                    directions,
                    args.levels,
                    rng,
                    args.state_probability,
                )
            )
            if not topology_edges:
                continue
            topology = FlowTopology(
                1 + max(max(edge[0], edge[1]) for edge in topology_edges),
                topology_edges,
                tuple(),
                tuple(),
            )
        elif args.random_cube:
            directions = tuple(rng.sample(pool, args.arms))
            topology_edges = tuple(
                random_cube_edges(
                    directions, rng, args.state_probability
                )
            )
            if not topology_edges:
                continue
            topology = FlowTopology(
                1 + max(max(edge[0], edge[1]) for edge in topology_edges),
                topology_edges,
                tuple(),
                tuple(),
            )
        elif args.punctured_cube:
            directions = tuple(rng.sample(pool, args.arms))
            topology_edges = tuple(punctured_cube_edges(directions))
            topology = FlowTopology(
                1 + max(max(edge[0], edge[1]) for edge in topology_edges),
                topology_edges,
                tuple(),
                tuple(),
            )
        elif args.hourglass:
            directions = (
                tuple(rng.sample(pool, rng.randint(2, 6))),
                tuple(rng.sample(pool, rng.randint(2, 6))),
            )
            topology_edges = tuple(hourglass(*directions))
            topology = FlowTopology(
                1 + max(max(edge[0], edge[1]) for edge in topology_edges),
                topology_edges,
                tuple(),
                tuple(),
            )
        else:
            directions = tuple(
                tuple(rng.sample(pool, 2)) for _ in range(args.arms)
            )
            topology = (
                capped_bowtie(directions, rng)
                if args.aggregate_cap
                else bowtie(directions, rng)
            )
        result = optimize(
            topology,
            args.minimum_bump,
            args.starts,
            args.iterations,
            rng,
            fixed_bump=args.fixed_bump,
            tangent_player=args.tangent_player,
        )
        if result is None:
            continue
        row = {
            "case": case,
            "directions": directions,
            "edges": [list(edge) for edge in topology.edges],
            "margin": result["margin"],
        }
        rows.append(row)
        if best is None or result["margin"] < best["margin"]:
            best = row | {"games": result["games"]}
            print(json.dumps({"best": row}), flush=True)
        if result["margin"] < -1e-8:
            break
    payload = {
        "status": (
            "negative_margin_found"
            if best is not None and best["margin"] < -1e-8
            else "none_found"
        ),
        "configuration": vars(args) | {"output": str(args.output)},
        "best": best,
        "rows": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
