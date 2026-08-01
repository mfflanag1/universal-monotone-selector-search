#!/usr/bin/env python3
"""Enumerate the complete binary 3-of-6 point-generated exact-game domain."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from array import array
from collections import defaultdict
from pathlib import Path

sys.path.append(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)

from n5_sparse_family_lp import sparse_family_slack_float


N = 6
GRAND = (1 << N) - 1
ALL_COALITIONS = (1 << (GRAND + 1)) - 1


class DisjointSet:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = bytearray(size)

    def find(self, item: int) -> int:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: int, right: int) -> None:
        left = self.find(left)
        right = self.find(right)
        if left == right:
            return
        if self.rank[left] < self.rank[right]:
            left, right = right, left
        self.parent[right] = left
        if self.rank[left] == self.rank[right]:
            self.rank[left] += 1


def point_planes(point: int) -> tuple[int, int, int]:
    planes = [0, 0, 0]
    for coalition in range(GRAND + 1):
        value = (point & coalition).bit_count()
        for threshold in range(1, 4):
            if value >= threshold:
                planes[threshold - 1] |= 1 << coalition
    return tuple(planes)  # type: ignore[return-value]


def game_from_planes(planes: tuple[int, int, int]) -> tuple[int, ...]:
    return tuple(
        sum((plane >> coalition) & 1 for plane in planes)
        for coalition in range(GRAND + 1)
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--solve-components", type=int, default=20)
    parser.add_argument("--max-solve-nodes", type=int, default=100000)
    args = parser.parse_args()

    points = [sum(1 << player for player in triple) for triple in itertools.combinations(range(N), 3)]
    point_data = [point_planes(point) for point in points]
    state_count = 1 << len(points)
    planes = [array("Q", [ALL_COALITIONS]), array("Q", [ALL_COALITIONS]), array("Q", [ALL_COALITIONS])]
    state_game = array("I", [0])
    game_index: dict[tuple[int, int, int], int] = {}
    game_planes: list[tuple[int, int, int]] = []

    for state in range(1, state_count):
        bit = state & -state
        point_index = bit.bit_length() - 1
        rest = state ^ bit
        current = tuple(
            planes[level][rest] & point_data[point_index][level]
            for level in range(3)
        )
        for level in range(3):
            planes[level].append(current[level])
        game = game_index.get(current)
        if game is None:
            game = len(game_planes)
            game_index[current] = game
            game_planes.append(current)
        state_game.append(game)
        if state % 100000 == 0:
            print(json.dumps({"phase": "states", "state": state, "games": len(game_planes)}), flush=True)

    dsu = DisjointSet(len(game_planes))
    edges: set[int] = set()
    edge_base = len(game_planes)
    for state in range(1, state_count):
        lower = state_game[state]
        remaining = state
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            successor = state ^ bit
            if not successor:
                continue
            upper = state_game[successor]
            if lower == upper:
                continue
            changed = (
                (planes[0][state] ^ planes[0][successor])
                | (planes[1][state] ^ planes[1][successor])
                | (planes[2][state] ^ planes[2][successor])
            )
            if changed.bit_count() != 1:
                continue
            coalition = changed.bit_length() - 1
            if coalition in (0, GRAND):
                continue
            encoded = (lower * edge_base + upper) * (GRAND + 1) + coalition
            if encoded not in edges:
                edges.add(encoded)
                dsu.union(lower, upper)
        if state % 100000 == 0:
            print(json.dumps({"phase": "edges", "state": state, "edges": len(edges)}), flush=True)

    components: dict[int, list[int]] = defaultdict(list)
    for game in range(len(game_planes)):
        components[dsu.find(game)].append(game)
    ordered = sorted(components.values(), key=len, reverse=True)
    edge_rows = []
    for encoded in edges:
        pair, coalition = divmod(encoded, GRAND + 1)
        lower, upper = divmod(pair, edge_base)
        edge_rows.append((lower, upper, coalition))
    component_edges: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for edge in edge_rows:
        component_edges[dsu.find(edge[0])].append(edge)

    summaries = []
    best = None
    for component_number, nodes in enumerate(ordered[: args.solve_components]):
        root = dsu.find(nodes[0])
        global_edges = component_edges[root]
        row = {
            "component": component_number,
            "nodes": len(nodes),
            "edges": len(global_edges),
            "cycle_rank": len(global_edges) - len(nodes) + 1,
            "margin": None,
        }
        if len(nodes) <= args.max_solve_nodes and global_edges:
            local = {node: index for index, node in enumerate(nodes)}
            games = [game_from_planes(game_planes[node]) for node in nodes]
            local_edges = [
                (local[lower], local[upper], coalition)
                for lower, upper, coalition in global_edges
            ]
            margin = sparse_family_slack_float(games, local_edges, N, "highs-ipm")
            row["margin"] = margin
            if margin is not None and (best is None or margin < best[0]):
                best = (margin, games, local_edges)
        summaries.append(row)
        print(json.dumps(row), flush=True)

    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "scope": "all nonempty subsets of the twenty binary 3-of-6 allocations; direct legal point deletions",
        "point_count": len(points),
        "state_count": state_count - 1,
        "distinct_game_count": len(game_planes),
        "legal_edge_count": len(edges),
        "component_count": len(ordered),
        "components": summaries,
        "best": None if best is None else {"margin": best[0], "nodes": len(best[1]), "edges": len(best[2])},
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if best is not None and best[0] < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
