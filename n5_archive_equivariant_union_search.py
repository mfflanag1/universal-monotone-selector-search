#!/usr/bin/env python3
"""Union every independently exact archived game and solve its symmetry quotient."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_equivariant_quotient_lp import quotient_family, solve_quotient
from n5_facet_search import load_facets


def add_induced_edges_fast(games, edges):
    grand = 31
    encoded = [
        tuple((value.numerator, value.denominator) for value in game)
        for game in games
    ]
    for coalition in range(1, grand):
        groups = {}
        for index, game in enumerate(encoded):
            signature = game[:coalition] + game[coalition + 1 :]
            groups.setdefault(signature, []).append(index)
        for indices in groups.values():
            if len(indices) < 2:
                continue
            ordered = sorted(indices, key=lambda index: games[index][coalition])
            for position, lower in enumerate(ordered):
                for upper in ordered[position + 1 :]:
                    edges.add((lower, upper, coalition))


def candidate_games(payload):
    family = payload.get("family")
    if isinstance(family, dict) and isinstance(family.get("games"), list):
        return family["games"]
    if isinstance(payload.get("games"), list):
        return payload["games"]
    return []


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facets, _ = load_facets()
    games = []
    game_set = set()
    accepted_files = 0
    rejected_games = 0
    for path in sorted(args.archive_root.glob("*.json")):
        try:
            payload = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        candidates = candidate_games(payload)
        accepted_here = False
        for raw_game in candidates:
            if not isinstance(raw_game, list):
                continue
            try:
                game = tuple(Fraction(str(value)) for value in raw_game)
            except (ValueError, ZeroDivisionError):
                continue
            if len(game) != 32 or game[0] != 0:
                continue
            accepted_here = True
            if game not in game_set:
                game_set.add(game)
                games.append(game)
        accepted_files += int(accepted_here)
    print(json.dumps({"phase": "collected", "distinct_games": len(games)}), flush=True)
    preliminary_edges = set()
    add_induced_edges_fast(games, preliminary_edges)
    preliminary_incident = {
        node for edge in preliminary_edges for node in edge[:2]
    }
    exact_games = []
    facet_matrix = np.asarray(facets, dtype=float)
    incident_list = sorted(preliminary_incident)
    for start in range(0, len(incident_list), 1000):
        indices = incident_list[start : start + 1000]
        values = np.asarray(
            [[float(value) for value in games[index]] for index in indices]
        )
        valid = np.min(values @ facet_matrix.T, axis=1) >= -1e-9
        for index, accepted in zip(indices, valid, strict=True):
            if accepted:
                exact_games.append(games[index])
            else:
                rejected_games += 1
    games = exact_games
    edges = set()
    add_induced_edges_fast(games, edges)
    incident = sorted({node for edge in edges for node in edge[:2]})
    reindex = {old: new for new, old in enumerate(incident)}
    active_games = [games[index] for index in incident]
    active_edges = {
        (reindex[lower], reindex[upper], coalition)
        for lower, upper, coalition in edges
    }
    print(
        json.dumps(
            {
                "accepted_files": accepted_files,
                "distinct_exact_games": len(games),
                "rejected_games": rejected_games,
                "active_games": len(active_games),
                "edges": len(active_edges),
            }
        ),
        flush=True,
    )
    margin = None
    orbit_games = 0
    quotient_edges = 0
    best_component = None
    if active_edges:
        canonical_games, quotient_edge_set, stabilizers = quotient_family(
            active_games, active_edges, 5
        )
        orbit_games = len(canonical_games)
        quotient_edges = len(quotient_edge_set)
        parent = list(range(orbit_games))

        def find(node):
            while parent[node] != node:
                parent[node] = parent[parent[node]]
                node = parent[node]
            return node

        def union(left, right):
            left = find(left)
            right = find(right)
            if left != right:
                parent[right] = left

        for lower, _lower_player, upper, _upper_player in quotient_edge_set:
            union(lower, upper)
        components = {}
        for node in range(orbit_games):
            components.setdefault(find(node), []).append(node)
        ordered_components = sorted(
            components.values(), key=len, reverse=True
        )
        for component_index, nodes in enumerate(ordered_components, start=1):
            node_set = set(nodes)
            local_index = {node: index for index, node in enumerate(nodes)}
            local_edges = {
                (
                    local_index[lower],
                    lower_player,
                    local_index[upper],
                    upper_player,
                )
                for lower, lower_player, upper, upper_player in quotient_edge_set
                if lower in node_set and upper in node_set
            }
            local_games = [canonical_games[node] for node in nodes]
            local_stabilizers = [stabilizers[node] for node in nodes]
            local_margin = solve_quotient(
                local_games, local_edges, local_stabilizers, 5
            )
            if margin is None or local_margin < margin:
                margin = local_margin
                best_component = (local_games, local_edges)
                print(
                    json.dumps(
                        {
                            "phase": "component",
                            "component": component_index,
                            "components": len(ordered_components),
                            "games": len(local_games),
                            "edges": len(local_edges),
                            "margin": margin,
                        }
                    ),
                    flush=True,
                )
            if local_margin < -1e-8:
                break
    found = margin is not None and margin < -1e-8
    payload = {
        "status": (
            "equivariant_obstruction_found"
            if found
            else "no_equivariant_obstruction_found"
        ),
        "archive_root": str(args.archive_root),
        "accepted_files": accepted_files,
        "distinct_exact_games": len(games),
        "rejected_games": rejected_games,
        "active_game_count": len(active_games),
        "edge_count": len(active_edges),
        "game_orbit_count": orbit_games,
        "quotient_edge_count": quotient_edges,
        "margin": margin,
        "quotient_component_games": (
            [list(map(str, game)) for game in best_component[0]]
            if found
            else None
        ),
        "quotient_component_edges": (
            [list(edge) for edge in sorted(best_component[1])]
            if found
            else None
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: value for key, value in payload.items() if key not in {"quotient_component_games", "quotient_component_edges"}}, indent=2))
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
