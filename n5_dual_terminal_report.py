#!/usr/bin/env python3
"""Report the terminal divergence topology of an exact margin dual."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import lcm
from pathlib import Path


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.archive.read_text())
    n = int(payload["n"])
    game_count = len(payload["family"]["games"])
    dual = payload["exact_margin_dual_certificate"]
    monotonicity_rows = [
        active
        for active in dual["active_inequalities"]
        if active["row"][0] == "monotonicity"
    ]
    scale = lcm(
        *(
            F(active["weight"]).denominator
            for active in monotonicity_rows
        )
    )
    divergence = [[F(0)] * n for _ in range(game_count)]
    active_edges = []
    player_flow: list[dict[tuple[int, int], F]] = [
        {} for _ in range(n)
    ]
    for active in monotonicity_rows:
        _, lower, upper, coalition, player = active["row"]
        flow = -F(active["weight"])
        divergence[int(lower)][int(player)] += flow
        divergence[int(upper)][int(player)] -= flow
        key = (int(lower), int(upper))
        player_flow[int(player)][key] = (
            player_flow[int(player)].get(key, F(0)) + flow * scale
        )
        active_edges.append(
            {
                "lower": int(lower),
                "upper": int(upper),
                "coalition": int(coalition),
                "player": int(player),
                "scaled_flow": str(flow * scale),
            }
        )

    terminals = []
    for node, vector in enumerate(divergence):
        if not any(vector):
            continue
        nonzero = {value for value in vector if value}
        indicator = len(nonzero) == 1
        coefficient = next(iter(nonzero)) if indicator else None
        block = (
            sum(
                1 << player
                for player, value in enumerate(vector)
                if value
            )
            if indicator
            else None
        )
        terminals.append(
            {
                "node": node,
                "scaled_divergence": [
                    str(value * scale) for value in vector
                ],
                "signed_indicator": indicator,
                "coefficient": (
                    None
                    if coefficient is None
                    else str(coefficient * scale)
                ),
                "coalition": block,
            }
        )

    paths = []
    for player in range(n):
        residual = dict(player_flow[player])
        supply = {
            node: vector[player] * scale
            for node, vector in enumerate(divergence)
            if vector[player] > 0
        }
        demand = {
            node: -vector[player] * scale
            for node, vector in enumerate(divergence)
            if vector[player] < 0
        }

        def find_path(source: int) -> list[int] | None:
            stack = [(source, [source])]
            while stack:
                node, path = stack.pop()
                if node in demand and demand[node] > 0:
                    return path
                for (lower, upper), flow in residual.items():
                    if lower == node and flow > 0 and upper not in path:
                        stack.append((upper, path + [upper]))
            return None

        for source in sorted(supply):
            while supply[source] > 0:
                path = find_path(source)
                if path is None:
                    raise RuntimeError(
                        f"cannot decompose player {player} flow at {source}"
                    )
                sink = path[-1]
                amount = min(supply[source], demand[sink])
                for lower, upper in zip(path, path[1:]):
                    amount = min(amount, residual[(lower, upper)])
                supply[source] -= amount
                demand[sink] -= amount
                for lower, upper in zip(path, path[1:]):
                    residual[(lower, upper)] -= amount
                paths.append(
                    {
                        "player": player,
                        "source": source,
                        "sink": sink,
                        "scaled_flow": str(amount),
                        "path": path,
                    }
                )
        if any(value for value in supply.values()) or any(
            value for value in demand.values()
        ):
            raise RuntimeError(f"unbalanced player {player} terminal flow")

    result = {
        "status": "exact_dual_terminal_report",
        "source": str(args.archive),
        "margin": payload["max_min_monotonicity_margin_exact"],
        "dual_objective": dual["objective"],
        "scale": scale,
        "game_count": game_count,
        "family_edge_count": len(payload["family"]["edges"]),
        "active_flow_row_count": len(active_edges),
        "active_core_row_count": sum(
            active["row"][0] == "core"
            for active in dual["active_inequalities"]
        ),
        "terminal_count": len(terminals),
        "all_signed_indicators": all(
            terminal["signed_indicator"] for terminal in terminals
        ),
        "terminals": terminals,
        "terminal_path_decomposition": paths,
        "active_flow_edges": active_edges,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
