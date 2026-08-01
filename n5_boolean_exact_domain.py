#!/usr/bin/env python3
"""Enumerate and solve the complete monotone Boolean five-player exact domain."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_facet_search import load_facets
from n5_point_deletion_family_search import connected_components
from n5_sparse_family_lp import sparse_family_slack_float


N = 5
GRAND = (1 << N) - 1


def antichains():
    comparable = []
    for coalition in range(1, GRAND + 1):
        mask = 0
        for other in range(1, GRAND + 1):
            if coalition & ~other == 0 or other & ~coalition == 0:
                mask |= 1 << (other - 1)
        comparable.append(mask)

    def visit(available: int, chosen: tuple[int, ...]):
        if not available:
            yield chosen
            return
        bit = available & -available
        index = bit.bit_length() - 1
        yield from visit(available ^ bit, chosen)
        yield from visit(available & ~comparable[index], chosen + (index + 1,))

    yield from visit((1 << GRAND) - 1, ())


def game_code(minimal_winning: tuple[int, ...]) -> int:
    code = 0
    for coalition in range(1, GRAND + 1):
        if any(generator & ~coalition == 0 for generator in minimal_winning):
            code |= 1 << coalition
    return code


def game_tuple(code: int) -> tuple[int, ...]:
    return tuple((code >> coalition) & 1 for coalition in range(GRAND + 1))


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facets, _types = load_facets()
    exact_codes = {0}
    monotone_count = 0
    for antichain in antichains():
        if not antichain:
            continue
        code = game_code(antichain)
        if not (code >> GRAND) & 1:
            continue
        monotone_count += 1
        game = game_tuple(code)
        if all(
            sum(facet[coalition] * game[coalition] for coalition in range(1, GRAND + 1)) >= 0
            for facet in facets
        ):
            exact_codes.add(code)
    codes = sorted(exact_codes)
    index = {code: position for position, code in enumerate(codes)}
    edges = set()
    for lower, code in enumerate(codes):
        for coalition in range(1, GRAND + 1):
            if code >> coalition & 1:
                continue
            upper_code = code | (1 << coalition)
            upper = index.get(upper_code)
            if upper is not None:
                edges.add((lower, upper, coalition))
    components = connected_components(len(codes), edges)
    rows = []
    best = None
    for component, (nodes, local_edges) in enumerate(components):
        games = [game_tuple(codes[node]) for node in nodes]
        margin = sparse_family_slack_float(games, local_edges, N, "highs-ipm")
        row = {
            "component": component,
            "nodes": len(nodes),
            "edges": len(local_edges),
            "cycle_rank": len(local_edges) - len(nodes) + 1,
            "margin": margin,
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
        if margin is not None and (best is None or margin < best[0]):
            best = (margin, games, local_edges)
    payload = {
        "status": (
            "incompatible_exact_family_found"
            if best is not None and best[0] < -1e-8
            else "no_incompatible_exact_family_found"
        ),
        "scope": "complete monotone Boolean five-player exact-game domain, including the zero game",
        "monotone_game_count": monotone_count,
        "exact_game_count": len(codes),
        "legal_edge_count": len(edges),
        "component_count": len(components),
        "components": rows,
        "best": None
        if best is None
        else {
            "margin": best[0],
            "games": [list(game) for game in best[1]] if best[0] < -1e-8 else None,
            "edges": [list(edge) for edge in best[2]] if best[0] < -1e-8 else None,
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if best is not None and best[0] < -1e-8 else 0


if __name__ == "__main__":
    raise SystemExit(main())
