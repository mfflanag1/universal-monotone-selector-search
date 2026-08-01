#!/usr/bin/env python3
"""Screen meet generators by their effect on a compact active component."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from n5_archive_equivariant_union_search import add_induced_edges_fast
from n5_equivariant_quotient_lp import quotient_family, solve_quotient
from n5_minkowski_family_product_search import load_games


def component_nodes(game_count, edges):
    parent = list(range(game_count))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for lower, _lower_player, upper, _upper_player in edges:
        left = find(lower)
        right = find(upper)
        if left != right:
            parent[right] = left
    groups = {}
    for node in range(game_count):
        groups.setdefault(find(node), []).append(node)
    return list(groups.values())


def local_problem(games, edges, stabilizers, nodes):
    node_set = set(nodes)
    reindex = {node: index for index, node in enumerate(nodes)}
    return (
        [games[node] for node in nodes],
        {
            (reindex[lower], lower_player, reindex[upper], upper_player)
            for lower, lower_player, upper, upper_player in edges
            if lower in node_set and upper in node_set
        },
        [stabilizers[node] for node in nodes],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base = load_games(args.base)
    candidates = load_games(args.candidates)
    if args.limit is not None:
        candidates = candidates[: args.limit]
    base_grand = {game[-1] for game in base}
    candidate_grand = {game[-1] for game in candidates}
    if len(base_grand) != 1 or len(candidate_grand) != 1:
        raise ValueError("each input family must have one grand worth")
    scale = next(iter(base_grand)) / next(iter(candidate_grand))
    candidates = [
        tuple(scale * value for value in game) for game in candidates
    ]

    base_edges = set()
    add_induced_edges_fast(base, base_edges)
    base_qgames, base_qedges, base_stabilizers = quotient_family(
        base, base_edges, 5
    )
    base_margin = solve_quotient(
        base_qgames, base_qedges, base_stabilizers, 5
    )
    base_qset = set(base_qgames)
    rows = []
    best_margin = base_margin
    best_family = None
    for candidate_index, candidate in enumerate(candidates):
        games = list(base)
        game_set = set(games)
        for game in base:
            meet = tuple(
                min(left_value, right_value)
                for left_value, right_value in zip(game, candidate)
            )
            if meet not in game_set:
                game_set.add(meet)
                games.append(meet)
        edges = set()
        add_induced_edges_fast(games, edges)
        qgames, qedges, stabilizers = quotient_family(games, edges, 5)
        parent_indices = {
            index for index, game in enumerate(qgames) if game in base_qset
        }
        affected = []
        for nodes in component_nodes(len(qgames), qedges):
            if parent_indices.intersection(nodes):
                affected.extend(nodes)
        local = local_problem(qgames, qedges, stabilizers, affected)
        unchanged = (
            len(local[0]) == len(base_qgames)
            and len(local[1]) == len(base_qedges)
        )
        margin = base_margin if unchanged else solve_quotient(*local, 5)
        row = {
            "candidate": candidate_index,
            "raw_games": len(games),
            "component_games": len(local[0]),
            "component_edges": len(local[1]),
            "margin": margin,
        }
        rows.append(row)
        if not unchanged or candidate_index % 25 == 0:
            print(json.dumps(row), flush=True)
        if margin < best_margin - 1e-12:
            best_margin = margin
            best_family = (games, edges, candidate_index)
            print(json.dumps({"phase": "best", **row}), flush=True)
        if margin < -1e-8:
            break

    output = {
        "status": (
            "equivariant_obstruction_found"
            if best_margin < -1e-8
            else "no_equivariant_obstruction_found"
        ),
        "base": str(args.base),
        "candidates": str(args.candidates),
        "base_margin": base_margin,
        "best_margin": best_margin,
        "best_candidate": None if best_family is None else best_family[2],
        "rows": rows,
    }
    if best_family is not None:
        output["family"] = {
            "games": [list(map(str, game)) for game in best_family[0]],
            "edges": [list(edge) for edge in sorted(best_family[1])],
        }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    return 1 if best_margin < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
