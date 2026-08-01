#!/usr/bin/env python3
"""Join exact families across grand-worth sections by all induced bump edges."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from family_search import Family, serial_family
from n5_sparse_dual_support import dual_envelope_support
from n5_sparse_family_lp import sparse_family_result


F = Fraction


def induced_edges(
    games: list[tuple[Fraction, ...]],
    grand: int,
) -> list[tuple[int, int, int]]:
    edges = set()
    for coalition in range(1, grand + 1):
        groups: dict[tuple[Fraction, ...], list[int]] = defaultdict(
            list
        )
        for index, game in enumerate(games):
            key = game[:coalition] + game[coalition + 1 :]
            groups[key].append(index)
        for indices in groups.values():
            ordered = sorted(
                indices, key=lambda index: games[index][coalition]
            )
            for position, lower in enumerate(ordered):
                for upper in ordered[position + 1 :]:
                    if (
                        games[lower][coalition]
                        < games[upper][coalition]
                    ):
                        edges.add((lower, upper, coalition))
    return sorted(edges)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--support-output", type=Path)
    parser.add_argument("--support-threshold", type=float, default=1e-9)
    parser.add_argument("--common-grand", type=F)
    parser.add_argument(
        "--method",
        choices=("highs", "highs-ds", "highs-ipm"),
        default="highs-ipm",
    )
    args = parser.parse_args()

    payloads = [json.loads(path.read_text()) for path in args.inputs]
    n = int(payloads[0]["n"])
    if any(int(payload["n"]) != n for payload in payloads):
        raise ValueError("player counts differ")
    grand = (1 << n) - 1
    raw_games = list(
        dict.fromkeys(
            tuple(F(value) for value in game)
            for payload in payloads
            for game in payload["family"]["games"]
        )
    )
    if args.common_grand is not None:
        if any(
            game[grand] > args.common_grand for game in raw_games
        ):
            raise ValueError("common grand is below a source grand worth")
        games = [
            game[:grand] + (args.common_grand,)
            for game in raw_games
        ]
        games = list(dict.fromkeys(games))
    else:
        games = raw_games
    edges = induced_edges(games, grand)
    print(
        json.dumps(
            {
                "stage": "induced_family_complete",
                "node_count": len(games),
                "edge_count": len(edges),
            }
        ),
        flush=True,
    )
    solve_result = sparse_family_result(
        games, edges, n, method=args.method
    )
    if not solve_result.success:
        raise RuntimeError(solve_result.message)
    margin = float(solve_result.x[len(games) * n])
    family = Family(
        games,
        edges,
        [0] * len(games),
        "n5_cross_grand_union",
    )
    aggregate_edges = sum(
        coalition == grand for _, _, coalition in edges
    )
    if args.support_output is not None:
        support = dual_envelope_support(
            games,
            edges,
            n,
            solve_result,
            str(args.output),
            args.support_threshold,
            args.method,
        )
        args.support_output.write_text(
            json.dumps(support, indent=2) + "\n"
        )
        print(
            json.dumps(
                {
                    "stage": "support_written",
                    "path": str(args.support_output),
                    "node_count": support["node_count"],
                    "edge_count": support["edge_count"],
                    "max_min_monotonicity_margin_float": support[
                        "max_min_monotonicity_margin_float"
                    ],
                }
            ),
            flush=True,
        )
    result = {
        "status": "n5_cross_grand_union_complete",
        "n": n,
        "sources": [str(path) for path in args.inputs],
        "node_count": len(games),
        "edge_count": len(edges),
        "aggregate_edge_count": aggregate_edges,
        "max_min_monotonicity_margin_float": margin,
        "family": serial_family(family, margin, None),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "node_count",
                    "edge_count",
                    "aggregate_edge_count",
                    "max_min_monotonicity_margin_float",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
