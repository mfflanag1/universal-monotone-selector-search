#!/usr/bin/env python3
"""Expand and dual-prune exact one-coordinate boundary neighbors."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(
    0,
    "/Users/maxf/projects/economics-research/game-theory/"
    "exact-game-monotone-selection/src",
)

from n5_equivariant_quotient_lp import quotient_family, solve_quotient
from n5_facet_search import load_facets
from n5_minkowski_family_product_search import load_games
from n5_prune_equivariant_dual_support import (
    components,
    dual_games,
    restrict_family,
)


def coordinate_limit(game, coalition, facets, increase):
    bounds = []
    for facet in facets:
        coefficient = facet[coalition]
        if (increase and coefficient >= 0) or (
            not increase and coefficient <= 0
        ):
            continue
        slack = sum(
            Fraction(value) * coefficient
            for value, coefficient in zip(game, facet, strict=True)
        )
        bound = (
            slack / -coefficient if increase else slack / coefficient
        )
        if bound < 0:
            raise ValueError("source game is outside the exact cone")
        bounds.append(bound)
    return min(bounds) if bounds else None


def best_component(games, edges, stabilizers):
    best = None
    for nodes in components(len(games), edges):
        local = restrict_family(games, edges, stabilizers, nodes)
        if not local[1]:
            continue
        margin = solve_quotient(*local, 5)
        if best is None or margin < best[0]:
            best = (margin, *local)
        if margin < -1e-8:
            break
    return best


def prune(games, edges, stabilizers, tolerance):
    rounds = []
    while True:
        margin, result = solve_quotient(
            games, edges, stabilizers, 5, return_result=True
        )
        active = dual_games(result, len(games), edges, tolerance)
        rounds.append(
            {
                "games": len(games),
                "edges": len(edges),
                "dual_games": len(active),
                "margin": margin,
            }
        )
        if not active or len(active) == len(games):
            return margin, games, edges, stabilizers, rounds
        candidate = restrict_family(games, edges, stabilizers, active)
        candidate_margin = solve_quotient(*candidate, 5)
        if abs(candidate_margin - margin) > 1e-8 * max(1.0, abs(margin)):
            rounds[-1]["rejected_candidate_margin"] = candidate_margin
            return margin, games, edges, stabilizers, rounds
        games, edges, stabilizers = candidate


def load_quotient(path):
    payload = json.loads(path.read_text())
    family = payload.get("family")
    if family is None:
        raise ValueError("input has no family")
    games = load_games(path)
    if family.get("quotient_edges") is not None:
        edges = {tuple(edge) for edge in family["quotient_edges"]}
        _canonical, _unused, stabilizers = quotient_family(games, set(), 5)
        if _canonical != games:
            raise ValueError("quotient input games are not canonical")
        return games, edges, stabilizers
    raw_edges = {tuple(edge) for edge in family.get("edges", [])}
    return quotient_family(games, raw_edges, 5)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    facets = load_facets()[0]
    games, edges, stabilizers = load_quotient(args.input)
    trace = []
    for round_index in range(1, args.rounds + 1):
        expanded = list(games)
        expanded_index = {game: index for index, game in enumerate(expanded)}
        boundary_edges = []
        for source, game in enumerate(games):
            for coalition in range(1, 31):
                for increase in (False, True):
                    delta = coordinate_limit(
                        game, coalition, facets, increase
                    )
                    if delta is None or delta <= 0:
                        continue
                    neighbor = list(game)
                    neighbor[coalition] += delta if increase else -delta
                    neighbor = tuple(neighbor)
                    target = expanded_index.get(neighbor)
                    if target is None:
                        target = len(expanded)
                        expanded_index[neighbor] = target
                        expanded.append(neighbor)
                    boundary_edges.append(
                        (source, target, coalition)
                        if increase
                        else (target, source, coalition)
                    )
        qgames, boundary_qedges, qstabilizers = quotient_family(
            expanded, boundary_edges, 5
        )
        if qgames[: len(games)] != games:
            raise RuntimeError("canonical source indices changed")
        qedges = set(boundary_qedges).union(edges)
        best = best_component(qgames, qedges, qstabilizers)
        if best is None:
            raise RuntimeError("boundary expansion has no component")
        margin, games, edges, stabilizers = best
        margin, games, edges, stabilizers, pruning = prune(
            games, edges, stabilizers, args.tolerance
        )
        row = {
            "round": round_index,
            "expanded_games": len(qgames),
            "expanded_edges": len(qedges),
            "margin": margin,
            "pruned_games": len(games),
            "pruned_edges": len(edges),
            "pruning": pruning,
        }
        trace.append(row)
        print(json.dumps(row), flush=True)
        args.output.write_text(
            json.dumps(
                {
                    "status": "boundary_closure_checkpoint",
                    "source": str(args.input),
                    "trace": trace,
                    "margin": margin,
                    "family": {
                        "games": [list(map(str, game)) for game in games],
                        "quotient_edges": [
                            list(edge) for edge in sorted(edges)
                        ],
                    },
                },
                indent=2,
            )
            + "\n"
        )
        if margin < -1e-8:
            break

    found = trace[-1]["margin"] < -1e-8
    output = {
        "status": (
            "equivariant_obstruction_found"
            if found
            else "no_equivariant_obstruction_found"
        ),
        "source": str(args.input),
        "trace": trace,
        "margin": trace[-1]["margin"],
        "family": {
            "games": [list(map(str, game)) for game in games],
            "quotient_edges": [list(edge) for edge in sorted(edges)],
        },
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
