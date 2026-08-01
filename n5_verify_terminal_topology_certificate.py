#!/usr/bin/env python3
"""Independently verify an exact terminal-topology relaxation certificate."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from n5_facet_search import load_facets


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    payload = json.loads(args.certificate.read_text())
    certificate = payload["exact_dual_certificate"]
    mixed_branches = payload.get("mixed_terminal_representations", {})
    if certificate is None:
        raise RuntimeError("archive has no exact dual certificate")
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    node_set = set(node_ids)
    node_index = {node: index for index, node in enumerate(node_ids)}
    variable_count = len(node_ids) * coalition_count

    def column(node: int, coalition: int) -> int:
        if node not in node_set or not 0 <= coalition <= grand:
            raise RuntimeError("certificate references an invalid variable")
        return node_index[node] * coalition_count + coalition

    protected_players: dict[tuple[int, int], set[int]] = defaultdict(set)
    for path in report["terminal_path_decomposition"]:
        protected_players[
            int(path["source"]), int(path["sink"])
        ].add(int(path["player"]))
    facets, _ = load_facets()
    stationarity = [F(0)] * variable_count
    dual_objective = F(0)

    for active in certificate["active_inequalities"]:
        metadata = active["row"]
        weight = F(active["weight"])
        if weight > 0:
            raise RuntimeError("an inequality dual weight is positive")
        row_type = metadata[0]
        if row_type == "exact_cone":
            _, node, facet_index = metadata
            facet = facets[int(facet_index)]
            for coalition, coefficient in enumerate(facet):
                if coefficient:
                    stationarity[column(int(node), coalition)] -= (
                        weight * F(coefficient)
                    )
        elif row_type == "game_monotonicity":
            _, node, coalition, successor = metadata
            node = int(node)
            coalition = int(coalition)
            successor = int(successor)
            difference = successor ^ coalition
            if (
                coalition & ~successor
                or difference == 0
                or difference & (difference - 1)
            ):
                raise RuntimeError("invalid game-monotonicity metadata")
            stationarity[column(node, coalition)] += weight
            stationarity[column(node, successor)] -= weight
        elif row_type == "directed_monotonicity":
            _, source, sink, coalition = metadata
            source = int(source)
            sink = int(sink)
            coalition = int(coalition)
            players = protected_players.get((source, sink))
            if not players:
                raise RuntimeError("directed row has no protected path")
            player_mask = sum(1 << player for player in players)
            if coalition & player_mask != player_mask:
                raise RuntimeError("directed row should be an invariance")
            stationarity[column(source, coalition)] += weight
            stationarity[column(sink, coalition)] -= weight
        else:
            raise RuntimeError(f"unknown inequality row type {row_type}")

    for active in certificate["active_equalities"]:
        metadata = active["row"]
        weight = F(active["weight"])
        row_type = metadata[0]
        if row_type == "empty":
            _, node = metadata
            stationarity[column(int(node), 0)] += weight
        elif row_type == "grand":
            _, node = metadata
            stationarity[column(int(node), grand)] += weight
            dual_objective += weight
        elif row_type == "protected_invariance":
            _, source, sink, coalition = metadata
            source = int(source)
            sink = int(sink)
            coalition = int(coalition)
            players = protected_players.get((source, sink))
            if not players:
                raise RuntimeError("invariance row has no protected path")
            player_mask = sum(1 << player for player in players)
            if coalition & player_mask == player_mask:
                raise RuntimeError("coalition is not protected by invariance")
            stationarity[column(source, coalition)] += weight
            stationarity[column(sink, coalition)] -= weight
        else:
            raise RuntimeError(f"unknown equality row type {row_type}")

    objective = [F(0)] * variable_count
    objective_constant = F(0)
    net_divergence = [F(0)] * n
    for terminal in terminals:
        node = int(terminal["node"])
        divergence = [F(value) for value in terminal["scaled_divergence"]]
        for player, value in enumerate(divergence):
            net_divergence[player] += value
        raw_coefficient = terminal.get("coefficient")
        if raw_coefficient is not None:
            coefficient = F(raw_coefficient)
            coalition = int(terminal["coalition"])
            if coefficient > 0:
                objective[column(node, coalition)] -= coefficient
            elif coefficient < 0:
                mass = -coefficient
                objective[column(node, grand ^ coalition)] -= mass
                objective_constant -= mass
        else:
            levels = sorted(set(divergence))
            if len(levels) == 2:
                low, high = levels
                coalition = sum(
                    1 << player
                    for player, value in enumerate(divergence)
                    if value == high
                )
                linear_terms = [(coalition, high - low), (grand, low)]
            else:
                selected = mixed_branches.get(str(node))
                if selected is None:
                    raise RuntimeError(
                        f"terminal {node} has no selected mixed branch"
                    )
                linear_terms = [
                    (grand, F(selected["mu"]))
                ] + [
                    (int(term["coalition"]), F(term["weight"]))
                    for term in selected["coalitions"]
                ]
                if any(weight < 0 for _, weight in linear_terms[1:]):
                    raise RuntimeError("mixed branch has a negative core weight")
                reconstructed = [F(selected["mu"])] * n
                for coalition, weight in linear_terms[1:]:
                    for player in range(n):
                        if coalition >> player & 1:
                            reconstructed[player] += weight
                if reconstructed != divergence:
                    raise RuntimeError(
                        f"terminal {node} branch does not represent divergence"
                    )
            for coalition, weight in linear_terms:
                objective[column(node, coalition)] -= weight
    if any(net_divergence):
        raise RuntimeError("terminal divergences do not conserve each player")
    if stationarity != objective:
        raise RuntimeError("exact dual stationarity failed")
    if dual_objective != objective_constant:
        raise RuntimeError("exact dual objective does not prove the zero bound")
    if F(certificate["variable_objective"]) != dual_objective:
        raise RuntimeError("reported variable objective is inconsistent")
    if F(certificate["lower_expectation_objective_with_constant"]) != 0:
        raise RuntimeError("reported lower-expectation bound is not zero")

    additive = [
        F(coalition.bit_count(), n) for coalition in range(coalition_count)
    ]
    if any(
        sum(
            (
                F(coefficient) * additive[coalition]
                for coalition, coefficient in enumerate(facet)
            ),
            F(0),
        )
        < 0
        for facet in facets
    ):
        raise RuntimeError("common additive primal violates the exact cone")
    primal_objective = F(0)
    for terminal in terminals:
        raw_coefficient = terminal.get("coefficient")
        if raw_coefficient is not None:
            coefficient = F(raw_coefficient)
            coalition = int(terminal["coalition"])
            if coefficient > 0:
                primal_objective += coefficient * additive[coalition]
            else:
                mass = -coefficient
                primal_objective += mass * (
                    additive[grand ^ coalition] - 1
                )
        else:
            divergence = [
                F(value) for value in terminal["scaled_divergence"]
            ]
            selected = mixed_branches.get(str(terminal["node"]))
            if len(set(divergence)) > 2 and selected is None:
                raise RuntimeError("missing mixed branch in primal check")
            primal_objective += sum(divergence, F(0)) * F(1, n)
    if primal_objective != 0:
        raise RuntimeError("common additive primal does not attain zero")

    print(
        "PASS: exact terminal-topology bound 0, "
        f"{len(node_ids)} terminals, {len(protected_players)} protected pairs, "
        f"{len(certificate['active_inequalities'])} inequality dual rows, "
        f"{len(certificate['active_equalities'])} equality dual rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
