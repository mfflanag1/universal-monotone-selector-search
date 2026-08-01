#!/usr/bin/env python3
"""Search large exact integer grids directly modulo all player permutations."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter, deque
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_boolean_exact_domain import GRAND, N
from n5_facet_search import load_facets
from n5_quaternary_reachable_domain import (
    exact_grand_two_games,
    local_monotone,
    reachable_section,
)


def permutation_tables():
    player_permutations = list(itertools.permutations(range(N)))
    inverse_players = []
    inverse_coalitions = []
    for permutation in player_permutations:
        inverse = [0] * N
        for old, new in enumerate(permutation):
            inverse[new] = old
        inverse_players.append(tuple(inverse))
        order = []
        for new_coalition in range(GRAND + 1):
            old_coalition = 0
            for new_player in range(N):
                if new_coalition >> new_player & 1:
                    old_coalition |= 1 << inverse[new_player]
            order.append(old_coalition)
        inverse_coalitions.append(order)
    return (
        player_permutations,
        inverse_players,
        np.asarray(inverse_coalitions, dtype=np.uint8),
    )


class OrbitDomain:
    def __init__(self, facet_matrix: np.ndarray, cap: int) -> None:
        self.facet_matrix = facet_matrix
        self.cap = cap
        (
            self.player_permutations,
            self.inverse_players,
            self.inverse_coalitions,
        ) = permutation_tables()
        self.alias: dict[bytes, tuple[int, tuple[int, ...]]] = {}
        self.games: list[bytes] = []
        self.player_variables: list[tuple[int, ...]] = []
        self.core_rows: dict[tuple[tuple[int, int], ...], int] = {}
        self.efficiency_rows: dict[tuple[tuple[int, int], ...], int] = {}
        self.next_variable = 0

    def orbit_rows(self, game: bytes) -> list[bytes]:
        array = np.frombuffer(game, dtype=np.uint8)
        matrix = array[self.inverse_coalitions]
        return [matrix[row].tobytes() for row in range(len(matrix))]

    def register(self, game: bytes) -> tuple[int, tuple[int, ...], bool]:
        existing = self.alias.get(game)
        if existing is not None:
            return existing[0], existing[1], False
        rows = self.orbit_rows(game)
        canonical = min(rows)
        canonicalizing = self.player_permutations[rows.index(canonical)]
        canonical_rows = self.orbit_rows(canonical)
        game_index = len(self.games)

        parent = list(range(N))

        def find(player: int) -> int:
            while parent[player] != player:
                parent[player] = parent[parent[player]]
                player = parent[player]
            return player

        def union(left: int, right: int) -> None:
            left = find(left)
            right = find(right)
            if left != right:
                parent[right] = left

        for permutation, row in zip(self.player_permutations, canonical_rows, strict=True):
            if row == canonical:
                for player in range(N):
                    union(player, permutation[player])
        root_variables: dict[int, int] = {}
        variables = []
        for player in range(N):
            root = find(player)
            if root not in root_variables:
                root_variables[root] = self.next_variable
                self.next_variable += 1
            variables.append(root_variables[root])
        self.games.append(canonical)
        self.player_variables.append(tuple(variables))

        for row, inverse in zip(canonical_rows, self.inverse_players, strict=True):
            self.alias.setdefault(row, (game_index, inverse))

        for coalition in range(1, GRAND):
            counts = Counter(
                variables[player]
                for player in range(N)
                if coalition >> player & 1
            )
            key = tuple(sorted(counts.items()))
            self.core_rows[key] = max(
                self.core_rows.get(key, canonical[coalition]), canonical[coalition]
            )
        counts = Counter(variables)
        key = tuple(sorted(counts.items()))
        previous = self.efficiency_rows.get(key)
        if previous is not None and previous != canonical[GRAND]:
            raise AssertionError("inconsistent quotient efficiency row")
        self.efficiency_rows[key] = canonical[GRAND]
        registered = self.alias[game]
        return registered[0], registered[1], True

    def exact(self, game: bytes) -> bool:
        values = np.frombuffer(game, dtype=np.uint8)[1:].astype(np.int16)
        return bool(np.all(self.facet_matrix @ values >= 0))

    def expand(self, seed_ids: list[int], max_orbits: int) -> bool:
        queue = deque(seed_ids)
        expanded = 0
        while queue:
            lower_id = queue.popleft()
            game = self.games[lower_id]
            score = self.facet_matrix @ np.frombuffer(game, dtype=np.uint8)[1:].astype(np.int16)
            for coalition in range(1, GRAND):
                current = game[coalition]
                for delta in (-1, 1):
                    value = current + delta
                    if (
                        value < 0
                        or value > self.cap
                        or not local_monotone(game, coalition, value)
                    ):
                        continue
                    if np.any(score + delta * self.facet_matrix[:, coalition - 1] < 0):
                        continue
                    candidate = bytearray(game)
                    candidate[coalition] = value
                    candidate_bytes = bytes(candidate)
                    existing = self.alias.get(candidate_bytes)
                    if existing is None and len(self.games) >= max_orbits:
                        return False
                    upper_id, _player_map, added = self.register(candidate_bytes)
                    if added:
                        queue.append(upper_id)
            expanded += 1
            if expanded % 1000 == 0:
                print(
                    json.dumps(
                        {
                            "phase": "orbit_bfs",
                            "expanded_orbits": expanded,
                            "orbits": len(self.games),
                            "labeled_games": len(self.alias),
                            "queue": len(queue),
                        }
                    ),
                    flush=True,
                )
        return True

    def monotonicity_rows(self) -> set[tuple[int, int]]:
        rows: set[tuple[int, int]] = set()
        for lower_id, game in enumerate(self.games):
            for coalition in range(1, GRAND + 1):
                if game[coalition] >= self.cap:
                    continue
                successor = bytearray(game)
                successor[coalition] += 1
                upper = self.alias.get(bytes(successor))
                if upper is None:
                    continue
                upper_id, player_map = upper
                for player in range(N):
                    if coalition >> player & 1:
                        rows.add(
                            (
                                self.player_variables[lower_id][player],
                                self.player_variables[upper_id][player_map[player]],
                            )
                        )
        return rows

    def game_components(self) -> list[list[int]]:
        parent = list(range(len(self.games)))

        def find(game: int) -> int:
            while parent[game] != game:
                parent[game] = parent[parent[game]]
                game = parent[game]
            return game

        def union(left: int, right: int) -> None:
            left = find(left)
            right = find(right)
            if left != right:
                parent[right] = left

        for lower_id, game in enumerate(self.games):
            for coalition in range(1, GRAND + 1):
                if game[coalition] >= self.cap:
                    continue
                successor = bytearray(game)
                successor[coalition] += 1
                upper = self.alias.get(bytes(successor))
                if upper is not None:
                    union(lower_id, upper[0])
        groups: dict[int, list[int]] = {}
        for game in range(len(self.games)):
            groups.setdefault(find(game), []).append(game)
        return sorted(groups.values(), key=len, reverse=True)

    def solve(self, feasibility_only: bool = False):
        monotonicity_rows = self.monotonicity_rows()
        slack = self.next_variable
        row_indices = []
        columns = []
        coefficients = []
        rhs = []
        row = 0
        for key, worth in self.core_rows.items():
            for variable, count in key:
                row_indices.append(row)
                columns.append(variable)
                coefficients.append(-float(count))
            rhs.append(-float(worth))
            row += 1
        for lower, upper in monotonicity_rows:
            row_indices.extend((row, row, row))
            columns.extend((lower, upper, slack))
            coefficients.extend((1.0, -1.0, 1.0))
            rhs.append(0.0)
            row += 1
        inequalities = coo_matrix(
            (coefficients, (row_indices, columns)),
            shape=(row, self.next_variable + 1),
        ).tocsr()
        eq_rows = []
        eq_columns = []
        eq_coefficients = []
        equality_rhs = []
        for eq_row, (key, grand_worth) in enumerate(self.efficiency_rows.items()):
            for variable, count in key:
                eq_rows.append(eq_row)
                eq_columns.append(variable)
                eq_coefficients.append(float(count))
            equality_rhs.append(float(grand_worth))
        equalities = coo_matrix(
            (eq_coefficients, (eq_rows, eq_columns)),
            shape=(len(self.efficiency_rows), self.next_variable + 1),
        ).tocsr()
        objective = np.zeros(self.next_variable + 1)
        if not feasibility_only:
            objective[slack] = -1.0
        bounds = [(None, None)] * (self.next_variable + 1)
        if feasibility_only:
            bounds[slack] = (0.0, 0.0)
        result = linprog(
            objective,
            A_ub=inequalities,
            b_ub=np.asarray(rhs),
            A_eq=equalities,
            b_eq=np.asarray(equality_rhs),
            bounds=bounds,
            method="highs-ipm",
        )
        return result, len(monotonicity_rows)

    def integer_certificate(self, denominator: int, seconds: float = 300.0):
        monotonicity_rows = self.monotonicity_rows()
        model = cp_model.CpModel()
        numerators = [
            model.new_int_var(0, denominator * self.cap, f"x_{index}")
            for index in range(self.next_variable)
        ]
        for key, worth in self.core_rows.items():
            model.add(
                sum(count * numerators[index] for index, count in key)
                >= denominator * worth
            )
        for key, grand_worth in self.efficiency_rows.items():
            model.add(
                sum(count * numerators[index] for index, count in key)
                == denominator * grand_worth
            )
        for lower, upper in monotonicity_rows:
            model.add(numerators[upper] - numerators[lower] >= 1)
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = seconds
        solver.parameters.num_search_workers = 8
        status = solver.solve(model)
        payload = {
            "status": solver.status_name(status),
            "denominator": denominator,
            "margin": f"1/{denominator}",
            "variable_count": self.next_variable,
            "core_row_count": len(self.core_rows),
            "efficiency_row_count": len(self.efficiency_rows),
            "monotonicity_row_count": len(monotonicity_rows),
            "wall_time": solver.wall_time,
            "verified": status in (cp_model.OPTIMAL, cp_model.FEASIBLE),
        }
        if payload["verified"]:
            values = [solver.value(variable) for variable in numerators]
            payload["allocation_numerators"] = values
            payload["verified"] = all(
                sum(count * values[index] for index, count in key) >= denominator * worth
                for key, worth in self.core_rows.items()
            ) and all(
                sum(count * values[index] for index, count in key) == denominator * grand_worth
                for key, grand_worth in self.efficiency_rows.items()
            ) and all(
                values[upper] - values[lower] >= 1
                for lower, upper in monotonicity_rows
            )
        return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cap", type=int, default=4)
    parser.add_argument("--max-orbits", type=int, default=20_000)
    parser.add_argument("--feasibility-only", action="store_true")
    parser.add_argument("--integer-denominator", type=int, default=0)
    parser.add_argument("--component-stats", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.cap != 4:
        raise ValueError("the current seed builder supports cap 4")

    facets, _types = load_facets()
    facet_matrix = np.asarray(facets, dtype=np.int16)[:, 1:]
    grand_two = exact_grand_two_games(facet_matrix)
    grand_three, complete = reachable_section(
        facet_matrix, grand_two, 3, 500_000
    )
    if not complete:
        raise AssertionError("grand-3 seed enumeration was unexpectedly truncated")

    domain = OrbitDomain(facet_matrix, args.cap)
    for game in itertools.chain(sorted(grand_two), grand_three):
        domain.register(game)
    print(
        json.dumps(
            {
                "phase": "lower_sections_registered",
                "orbits": len(domain.games),
                "labeled_games": len(domain.alias),
            }
        ),
        flush=True,
    )
    seed_ids = []
    for game in grand_three:
        raised = bytearray(game)
        raised[GRAND] = args.cap
        game_id, _player_map, added = domain.register(bytes(raised))
        if added:
            seed_ids.append(game_id)
    print(
        json.dumps(
            {
                "phase": "cap4_seeds_registered",
                "orbits": len(domain.games),
                "labeled_games": len(domain.alias),
                "seed_orbits": len(seed_ids),
            }
        ),
        flush=True,
    )
    del grand_three
    complete = domain.expand(seed_ids, args.max_orbits)
    if args.component_stats:
        components = domain.game_components()
        rows = [
            {
                "component": index,
                "game_orbits": len(component),
                "allocation_variable_orbits": len(
                    {
                        variable
                        for game in component
                        for variable in domain.player_variables[game]
                    }
                ),
            }
            for index, component in enumerate(components)
        ]
        payload = {
            "status": "component_stats_complete",
            "bfs_complete": complete,
            "orbit_count": len(domain.games),
            "labeled_game_count": len(domain.alias),
            "component_count": len(components),
            "components": rows,
        }
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
        print(json.dumps(payload, indent=2))
        return 0
    if args.integer_denominator:
        certificate = domain.integer_certificate(args.integer_denominator)
        payload = {
            "status": (
                "exact_positive_certificate_found"
                if certificate["verified"]
                else "integer_certificate_not_found"
            ),
            "scope": "permutation-orbit closure of the reachable exact integer grid through grand worth 4",
            "bfs_complete": complete,
            "orbit_count": len(domain.games),
            "labeled_game_count": len(domain.alias),
            "allocation_variable_orbits": domain.next_variable,
            "core_row_orbits": len(domain.core_rows),
            "efficiency_row_orbits": len(domain.efficiency_rows),
            "integer_primal_certificate": certificate,
        }
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
        print(
            json.dumps(
                {key: value for key, value in payload.items() if key != "integer_primal_certificate"}
                | {"certificate": {key: value for key, value in certificate.items() if key != "allocation_numerators"}},
                indent=2,
            )
        )
        return 0 if certificate["verified"] else 2
    result, monotonicity_row_count = domain.solve(args.feasibility_only)
    if not result.success:
        if args.feasibility_only:
            payload = {
                "status": "weak_incompatibility_detected_float",
                "scope": "permutation-orbit closure of the reachable exact integer grid through grand worth 4",
                "bfs_complete": complete,
                "orbit_count": len(domain.games),
                "labeled_game_count": len(domain.alias),
                "message": result.message,
            }
            args.output.write_text(json.dumps(payload, indent=2) + "\n")
            print(json.dumps(payload, indent=2))
            return 1
        raise RuntimeError(result.message)
    margin = float(result.x[-1])
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if margin < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "scope": "permutation-orbit closure of the reachable exact integer grid through grand worth 4",
        "bfs_complete": complete,
        "orbit_count": len(domain.games),
        "labeled_game_count": len(domain.alias),
        "allocation_variable_orbits": domain.next_variable,
        "core_row_orbits": len(domain.core_rows),
        "efficiency_row_orbits": len(domain.efficiency_rows),
        "monotonicity_row_orbits": monotonicity_row_count,
        "margin": margin,
        "feasibility_only": args.feasibility_only,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if margin < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
