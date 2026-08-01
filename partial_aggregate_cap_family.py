#!/usr/bin/env python3
"""Attach selected aggregate predecessors to a proper-bump exact family."""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

sys.path.append(
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src"
)

from n5_sparse_family_lp import sparse_family_slack_float


def build(games, edges, n, delta, capped):
    grand = (1 << n) - 1
    result_games = []
    for game in games:
        upper = list(game)
        upper[grand] += delta
        result_games.append(upper)
    result_edges = [tuple(edge) for edge in edges]
    for node in capped:
        lower = list(games[node])
        lower_index = len(result_games)
        result_games.append(lower)
        result_edges.append((lower_index, node, grand))
    return result_games, result_edges


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--grand-delta", type=float, default=1.0)
    parser.add_argument("--trials", type=int, default=20)
    parser.add_argument("--cap-fraction", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=20260731)
    args = parser.parse_args()
    source = json.loads(args.input.read_text())
    games = source["best"]["games"]
    edges = source["best"]["edges"]
    n = len(source["pattern"])
    rng = random.Random(args.seed)
    rows = []
    best = None
    node_count = len(games)
    for trial in range(args.trials):
        if trial == 0 or args.cap_fraction >= 1:
            capped = list(range(node_count))
        else:
            count = max(1, round(args.cap_fraction * node_count))
            capped = sorted(rng.sample(range(node_count), count))
        candidate_games, candidate_edges = build(
            games, edges, n, args.grand_delta, capped
        )
        margin = sparse_family_slack_float(
            candidate_games, candidate_edges, n, "highs-ipm"
        )
        row = {"trial": trial, "caps": len(capped), "margin": margin}
        rows.append(row)
        print(json.dumps(row), flush=True)
        if margin is not None and (best is None or margin < best[0]):
            best = (margin, capped, candidate_games, candidate_edges)
        if margin is not None and margin < -1e-8:
            break
        if args.cap_fraction >= 1:
            break
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "source": str(args.input),
        "n": n,
        "grand_delta": args.grand_delta,
        "rows": rows,
        "best": None
        if best is None
        else {
            "margin": best[0],
            "capped_nodes": best[1],
            "games": best[2] if best[0] < -1e-8 else None,
            "edges": [list(edge) for edge in best[3]] if best[0] < -1e-8 else None,
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"status": payload["status"], "best_margin": None if best is None else best[0]}, indent=2))
    return 1 if best is not None and best[0] < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
