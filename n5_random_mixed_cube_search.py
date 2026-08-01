#!/usr/bin/env python3
"""Adversarial search on random non-product subcomplexes of an exact-game cube."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from n5_facet_search import search


def connected_states(dimension: int, count: int, rng: random.Random) -> set[int]:
    states = {0}
    frontier = [0]
    while len(states) < count:
        parent = rng.choice(frontier)
        missing = [bit for bit in range(dimension) if not parent >> bit & 1]
        if not missing:
            frontier.remove(parent)
            if not frontier:
                break
            continue
        child = parent | (1 << rng.choice(missing))
        states.add(child)
        frontier.append(child)
        if rng.random() < 0.35:
            for bit in range(dimension):
                neighbor = child ^ (1 << bit)
                if neighbor in states and neighbor not in frontier:
                    frontier.append(neighbor)
    return states


def induced_edges(
    states: set[int], directions: list[int]
) -> list[tuple[int, int, int]]:
    ordered = sorted(states)
    index = {state: node for node, state in enumerate(ordered)}
    edges = []
    for state in ordered:
        for bit, coalition in enumerate(directions):
            child = state | (1 << bit)
            if child != state and child in states:
                edges.append((index[state], index[child], coalition))
    return edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=100)
    parser.add_argument("--dimension", type=int, default=9)
    parser.add_argument("--nodes", type=int, default=64)
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=12)
    parser.add_argument("--bump", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    rows = []
    best = None
    for case in range(args.cases):
        proper = rng.sample(range(1, 31), args.dimension - 1)
        directions = proper + [31]
        rng.shuffle(directions)
        states = connected_states(
            args.dimension, min(args.nodes, 1 << args.dimension), rng
        )
        edges = induced_edges(states, directions)
        if not edges or 31 not in {edge[2] for edge in edges}:
            continue
        bumps_by_direction = {
            coalition: args.bump * rng.randint(1, 8) / 4
            for coalition in directions
        }
        bumps = [bumps_by_direction[coalition] for _, _, coalition in edges]
        try:
            record = search(
                edges,
                bumps,
                args.starts,
                args.iterations,
                args.seed + case,
                quiet=True,
            )
        except RuntimeError:
            continue
        margin = float(record["best_margin_float"])
        row = {
            "case": case,
            "node_count": len(states),
            "edge_count": len(edges),
            "directions": directions,
            "margin": margin,
        }
        rows.append(row)
        if best is None or margin < float(best["best_margin_float"]):
            best = record
            print(json.dumps({"best": row}), flush=True)
        elif case % 10 == 0:
            print(json.dumps(row), flush=True)
        if margin < -1e-8:
            break
    payload = {
        "status": "negative_found" if best and best["best_margin_float"] < -1e-8 else "no_negative_found",
        "configuration": vars(args) | {"output": str(args.output)},
        "rows": rows,
        "best": best,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
