#!/usr/bin/env python3
"""Independently verify an exact reduced perspective dual certificate."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any

from n5_facet_search import load_facets
from n5_mixed_terminal_branch_search import extreme_representations


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    certificate = json.loads(args.certificate.read_text())
    if certificate["status"] != "mixed_terminal_perspective_reduced_exact":
        raise RuntimeError("unexpected certificate status")

    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    node_ids = sorted(int(terminal["node"]) for terminal in terminals)
    node_set = set(node_ids)
    node_index = {node: index for index, node in enumerate(node_ids)}
    base_variable_count = len(node_ids) * coalition_count

    mixed_terminals = [
        terminal
        for terminal in terminals
        if terminal.get("coefficient") is None
        and len(set(terminal["scaled_divergence"])) > 2
    ]
    mixed_nodes = sorted(int(terminal["node"]) for terminal in mixed_terminals)
    representations = {
        int(terminal["node"]): extreme_representations(
            [F(value) for value in terminal["scaled_divergence"]], n
        )
        for terminal in mixed_terminals
    }

    choice_offsets: dict[int, int] = {}
    copy_offsets: dict[int, int] = {}
    next_variable = base_variable_count
    for terminal in mixed_terminals:
        node = int(terminal["node"])
        choice_offsets[node] = next_variable
        next_variable += len(representations[node])
        copy_offsets[node] = next_variable
        next_variable += len(representations[node]) * coalition_count
    variable_count = next_variable

    def base_column(node: int, coalition: int) -> int:
        if node not in node_set or not 0 <= coalition <= grand:
            raise RuntimeError("invalid base-game coordinate")
        return node_index[node] * coalition_count + coalition

    def choice_column(node: int, choice: int) -> int:
        if node not in representations or not 0 <= choice < len(
            representations[node]
        ):
            raise RuntimeError("invalid perspective choice coordinate")
        return choice_offsets[node] + choice

    def copy_column(node: int, choice: int, coalition: int) -> int:
        if not 0 <= coalition <= grand:
            raise RuntimeError("invalid perspective-copy coalition")
        choice_column(node, choice)
        return copy_offsets[node] + choice * coalition_count + coalition

    facets, _ = load_facets()
    protected_players: dict[tuple[int, int], set[int]] = defaultdict(set)
    for path in report["terminal_path_decomposition"]:
        protected_players[
            int(path["source"]), int(path["sink"])
        ].add(int(path["player"]))

    objective = [F(0)] * variable_count
    objective_constant = F(0)
    for terminal in terminals:
        node = int(terminal["node"])
        if node in representations:
            continue
        raw_coefficient = terminal.get("coefficient")
        if raw_coefficient is not None:
            coefficient = F(raw_coefficient)
            coalition = int(terminal["coalition"])
            if coefficient > 0:
                objective[base_column(node, coalition)] -= coefficient
            elif coefficient < 0:
                mass = -coefficient
                objective[base_column(node, grand ^ coalition)] -= mass
                objective_constant -= mass
        else:
            divergence = [F(value) for value in terminal["scaled_divergence"]]
            low, high = sorted(set(divergence))
            coalition = sum(
                1 << player
                for player, value in enumerate(divergence)
                if value == high
            )
            objective[base_column(node, coalition)] -= high - low
            objective[base_column(node, grand)] -= low
    for node, choices in representations.items():
        for choice, (mu, coalitions) in enumerate(choices):
            objective[choice_column(node, choice)] -= mu
            for coalition, weight in coalitions:
                objective[copy_column(node, choice, coalition)] -= weight

    if F(certificate["objective_constant"]) != objective_constant:
        raise RuntimeError("reported objective constant is incorrect")

    stationarity = [F(0)] * variable_count
    dual_objective = F(0)
    inequality_count = 0
    equality_count = 0

    def add(entries: dict[int, F], weight: F) -> None:
        for column, coefficient in entries.items():
            stationarity[column] += weight * coefficient

    seen_global_rows = set()
    for active in certificate["global_certificate"]["active_inequalities"]:
        metadata = tuple(active["row"])
        if metadata in seen_global_rows:
            raise RuntimeError("duplicate global dual row")
        seen_global_rows.add(metadata)
        weight = F(active["weight"])
        if weight >= 0:
            raise RuntimeError("global inequality dual is not negative")
        row_type = metadata[0]
        if row_type == "exact_cone":
            _, node, facet_index = metadata
            facet = facets[int(facet_index)]
            add(
                {
                    base_column(int(node), coalition): -F(coefficient)
                    for coalition, coefficient in enumerate(facet)
                    if coefficient
                },
                weight,
            )
        elif row_type == "game_monotonicity":
            _, node, coalition, successor = metadata
            coalition = int(coalition)
            successor = int(successor)
            if coalition & successor != coalition or (
                successor ^ coalition
            ).bit_count() != 1:
                raise RuntimeError("invalid game-monotonicity row")
            add(
                {
                    base_column(int(node), coalition): F(1),
                    base_column(int(node), successor): F(-1),
                },
                weight,
            )
        elif row_type == "directed_monotonicity":
            _, source, sink, coalition = metadata
            source = int(source)
            sink = int(sink)
            coalition = int(coalition)
            players = protected_players[source, sink]
            mask = sum(1 << player for player in players)
            if not players or coalition & mask != mask:
                raise RuntimeError("invalid directed-monotonicity row")
            add(
                {
                    base_column(source, coalition): F(1),
                    base_column(sink, coalition): F(-1),
                },
                weight,
            )
        else:
            raise RuntimeError(f"unsupported global inequality {row_type}")
        inequality_count += 1

    for active in certificate["global_certificate"]["active_equalities"]:
        metadata = tuple(active["row"])
        if metadata in seen_global_rows:
            raise RuntimeError("duplicate global dual row")
        seen_global_rows.add(metadata)
        weight = F(active["weight"])
        row_type = metadata[0]
        if row_type == "empty":
            _, node = metadata
            add({base_column(int(node), 0): F(1)}, weight)
        elif row_type == "grand":
            _, node = metadata
            add({base_column(int(node), grand): F(1)}, weight)
            dual_objective += weight
        elif row_type == "protected_invariance":
            _, source, sink, coalition = metadata
            source = int(source)
            sink = int(sink)
            coalition = int(coalition)
            players = protected_players[source, sink]
            mask = sum(1 << player for player in players)
            if not players or coalition & mask == mask:
                raise RuntimeError("invalid protected-invariance row")
            add(
                {
                    base_column(source, coalition): F(1),
                    base_column(sink, coalition): F(-1),
                },
                weight,
            )
        else:
            raise RuntimeError(f"unsupported global equality {row_type}")
        equality_count += 1

    anchors = certificate["anchors"]
    simplex = certificate["simplex"]
    if set(anchors) != {str(node) for node in mixed_nodes}:
        raise RuntimeError("anchor nodes are incomplete")
    if set(simplex) != {str(node) for node in mixed_nodes}:
        raise RuntimeError("simplex nodes are incomplete")
    for node in mixed_nodes:
        link_weights = [F(value) for value in anchors[str(node)]]
        if len(link_weights) != coalition_count:
            raise RuntimeError("incorrect perspective-link count")
        for coalition, weight in enumerate(link_weights):
            row = {base_column(node, coalition): F(1)}
            for choice in range(len(representations[node])):
                row[copy_column(node, choice, coalition)] = F(-1)
            add(row, weight)
            equality_count += bool(weight)
        simplex_weight = F(simplex[str(node)])
        add(
            {
                choice_column(node, choice): F(1)
                for choice in range(len(representations[node]))
            },
            simplex_weight,
        )
        dual_objective += simplex_weight
        equality_count += bool(simplex_weight)

    local_certificates = certificate["local_certificates"]
    if set(local_certificates) != {str(node) for node in mixed_nodes}:
        raise RuntimeError("local certificate nodes are incomplete")
    for node in mixed_nodes:
        blocks = local_certificates[str(node)]
        if len(blocks) != len(representations[node]):
            raise RuntimeError("local certificate choices are incomplete")
        for expected_choice, block in enumerate(blocks):
            if int(block["choice"]) != expected_choice:
                raise RuntimeError("local certificate choices are unordered")
            mu, coalitions = representations[node][expected_choice]
            reported_coalitions = [
                (int(entry["coalition"]), F(entry["weight"]))
                for entry in block["coalitions"]
            ]
            if F(block["mu"]) != mu or tuple(reported_coalitions) != coalitions:
                raise RuntimeError("mixed-terminal representation mismatch")
            seen_local_rows = set()
            for active in block["active_rows"]:
                metadata = tuple(active["row"])
                if metadata in seen_local_rows:
                    raise RuntimeError("duplicate local dual row")
                seen_local_rows.add(metadata)
                weight = F(active["weight"])
                row_type = metadata[0]
                if row_type == "perspective_copy_exact_cone":
                    _, row_node, choice, facet_index = metadata
                    if int(row_node) != node or int(choice) != expected_choice:
                        raise RuntimeError("local exact-cone row is misassigned")
                    if weight >= 0:
                        raise RuntimeError("local inequality dual is not negative")
                    facet = facets[int(facet_index)]
                    add(
                        {
                            copy_column(
                                node, expected_choice, coalition
                            ): -F(coefficient)
                            for coalition, coefficient in enumerate(facet)
                            if coefficient
                        },
                        weight,
                    )
                    inequality_count += 1
                elif row_type == "perspective_copy_game_monotonicity":
                    _, row_node, choice, coalition, successor = metadata
                    coalition = int(coalition)
                    successor = int(successor)
                    if (
                        int(row_node) != node
                        or int(choice) != expected_choice
                        or coalition & successor != coalition
                        or (successor ^ coalition).bit_count() != 1
                    ):
                        raise RuntimeError("invalid local monotonicity row")
                    if weight >= 0:
                        raise RuntimeError("local inequality dual is not negative")
                    add(
                        {
                            copy_column(
                                node, expected_choice, coalition
                            ): F(1),
                            copy_column(
                                node, expected_choice, successor
                            ): F(-1),
                        },
                        weight,
                    )
                    inequality_count += 1
                elif row_type == "perspective_choice_nonnegative":
                    _, row_node, choice = metadata
                    if int(row_node) != node or int(choice) != expected_choice:
                        raise RuntimeError("invalid choice-nonnegative row")
                    if weight >= 0:
                        raise RuntimeError("local inequality dual is not negative")
                    add(
                        {choice_column(node, expected_choice): F(-1)},
                        weight,
                    )
                    inequality_count += 1
                elif row_type == "perspective_copy_empty":
                    _, row_node, choice = metadata
                    if int(row_node) != node or int(choice) != expected_choice:
                        raise RuntimeError("invalid copy-empty row")
                    add(
                        {copy_column(node, expected_choice, 0): F(1)},
                        weight,
                    )
                    equality_count += 1
                elif row_type == "perspective_copy_grand":
                    _, row_node, choice = metadata
                    if int(row_node) != node or int(choice) != expected_choice:
                        raise RuntimeError("invalid copy-grand row")
                    add(
                        {
                            copy_column(
                                node, expected_choice, grand
                            ): F(1),
                            choice_column(node, expected_choice): F(-1),
                        },
                        weight,
                    )
                    equality_count += 1
                else:
                    raise RuntimeError(f"unsupported local row {row_type}")

    if stationarity != objective:
        residuals = [
            (index, expected - actual)
            for index, (expected, actual) in enumerate(
                zip(objective, stationarity, strict=True)
            )
            if expected != actual
        ]
        raise RuntimeError(
            f"exact stationarity failed in {len(residuals)} coordinates; "
            f"first={residuals[:3]}"
        )
    if dual_objective != objective_constant:
        raise RuntimeError("exact dual objective does not prove the zero bound")

    additive = [F(coalition.bit_count(), n) for coalition in range(coalition_count)]
    if additive[0] != 0 or additive[grand] != 1:
        raise RuntimeError("additive normalization failed")
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
        raise RuntimeError("additive game violates the exact cone")
    primal = [F(0)] * variable_count
    for node in node_ids:
        for coalition, value in enumerate(additive):
            primal[base_column(node, coalition)] = value
    for node in mixed_nodes:
        primal[choice_column(node, 0)] = F(1)
        for coalition, value in enumerate(additive):
            primal[copy_column(node, 0, coalition)] = value
    primal_objective = sum(
        (
            coefficient * value
            for coefficient, value in zip(objective, primal, strict=True)
        ),
        F(0),
    )
    if primal_objective != objective_constant:
        raise RuntimeError("common additive primal does not attain the bound")

    if inequality_count != (
        len(certificate["global_certificate"]["active_inequalities"])
        + int(certificate["local_inequality_row_count"])
    ):
        raise RuntimeError("reported inequality-row count is inconsistent")
    print(
        "PASS: exact m=8 perspective bound 0, "
        f"{variable_count} variables, {inequality_count} inequality rows, "
        f"{equality_count} nonzero equality rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
