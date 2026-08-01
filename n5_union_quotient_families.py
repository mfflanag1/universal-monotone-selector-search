#!/usr/bin/env python3
"""Union compact equivariant exact families and dual-prune the result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from n5_equivariant_boundary_closure import best_component, load_quotient, prune
from n5_equivariant_quotient_lp import quotient_family


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    games = []
    game_index = {}
    edges = set()
    rows = []
    for path in args.inputs:
        local_games, local_edges, _local_stabilizers = load_quotient(path)
        remap = []
        for game in local_games:
            target = game_index.get(game)
            if target is None:
                target = len(games)
                game_index[game] = target
                games.append(game)
            remap.append(target)
        for lower, lower_player, upper, upper_player in local_edges:
            edges.add(
                (
                    remap[lower],
                    lower_player,
                    remap[upper],
                    upper_player,
                )
            )
        rows.append(
            {
                "input": str(path),
                "games": len(local_games),
                "edges": len(local_edges),
            }
        )
    canonical, _unused, stabilizers = quotient_family(games, set(), 5)
    if canonical != games:
        raise RuntimeError("input quotient games are not canonical")
    best = best_component(games, edges, stabilizers)
    if best is None:
        raise RuntimeError("union has no comparison component")
    margin, games, edges, stabilizers = best
    margin, games, edges, stabilizers, pruning = prune(
        games, edges, stabilizers, args.tolerance
    )
    found = margin < -1e-8
    output = {
        "status": (
            "equivariant_obstruction_found"
            if found
            else "no_equivariant_obstruction_found"
        ),
        "inputs": rows,
        "union_games": len(game_index),
        "union_edges": sum(row["edges"] for row in rows),
        "margin": margin,
        "pruned_games": len(games),
        "pruned_edges": len(edges),
        "pruning": pruning,
        "family": {
            "games": [list(map(str, game)) for game in games],
            "quotient_edges": [list(edge) for edge in sorted(edges)],
        },
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: value for key, value in output.items() if key != "family"}, indent=2))
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
