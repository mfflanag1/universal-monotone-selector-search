#!/usr/bin/env python3
"""Solve a serialized family with the memory-safe sparse LP."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from n5_sparse_dual_support import dual_envelope_support
from n5_sparse_family_lp import sparse_family_result


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--method",
        choices=("highs", "highs-ds", "highs-ipm"),
        default="highs",
    )
    parser.add_argument("--disp", action="store_true")
    parser.add_argument("--no-crossover", action="store_true")
    parser.add_argument("--support-output", type=Path)
    parser.add_argument("--support-threshold", type=float, default=1e-9)
    parser.add_argument(
        "--dual-edge-weight-strategy",
        choices=("dantzig", "devex", "steepest", "steepest-devex"),
    )
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    n = int(payload["n"])
    raw = payload["family"]
    games = [
        tuple(F(value) for value in game)
        for game in raw["games"]
    ]
    edges = [
        (
            int(edge["lower"]),
            int(edge["upper"]),
            int(edge["coalition"]),
        )
        for edge in raw["edges"]
    ]
    solve_result = sparse_family_result(
        games,
        edges,
        n,
        method=args.method,
        dual_edge_weight_strategy=args.dual_edge_weight_strategy,
        disp=args.disp,
        run_crossover=False if args.no_crossover else None,
    )
    if not solve_result.success:
        raise RuntimeError(solve_result.message)
    margin = float(solve_result.x[len(games) * n])
    support = None
    if args.support_output is not None:
        support = dual_envelope_support(
            games,
            edges,
            n,
            solve_result,
            str(args.input),
            args.support_threshold,
            args.method,
        )
        args.support_output.write_text(json.dumps(support, indent=2) + "\n")
    result = {
        "status": "sparse_archive_solve_complete",
        "n": n,
        "source": str(args.input),
        "method": args.method,
        "dual_edge_weight_strategy": args.dual_edge_weight_strategy,
        "run_crossover": not args.no_crossover,
        "node_count": len(games),
        "edge_count": len(edges),
        "max_min_monotonicity_margin_float": margin,
        "support_output": (
            str(args.support_output) if args.support_output is not None else None
        ),
        "support_node_count": support["node_count"] if support is not None else None,
        "support_edge_count": support["edge_count"] if support is not None else None,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
