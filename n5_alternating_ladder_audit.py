#!/usr/bin/env python3
"""Recognize and exactly verify a complementary alternating-ladder dual."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any


F = Fraction


def coalition_members(coalition: int, n: int) -> tuple[int, ...]:
    return tuple(player for player in range(n) if coalition >> player & 1)


def recognize_orientation(
    sources: list[int],
    sinks: list[dict[str, Any]],
    assignments: dict[tuple[int, int], int],
    n: int,
    block_a: int,
    block_b: int,
) -> dict[str, Any] | None:
    grand = (1 << n) - 1
    grand_sinks = [sink for sink in sinks if sink["coalition"] == grand]
    sink_a = next(
        (sink for sink in sinks if sink["coalition"] == block_a), None
    )
    sink_b = next(
        (sink for sink in sinks if sink["coalition"] == block_b), None
    )
    if sink_a is None or sink_b is None:
        return None

    relations: list[dict[str, int]] = []
    for sink in grand_sinks:
        node = int(sink["node"])
        a_sources = {
            source
            for source in sources
            for player in coalition_members(block_a, n)
            if assignments.get((source, player)) == node
        }
        b_sources = {
            source
            for source in sources
            for player in coalition_members(block_b, n)
            if assignments.get((source, player)) == node
        }
        if len(a_sources) != 1 or len(b_sources) != 1:
            return None
        source_a = next(iter(a_sources))
        source_b = next(iter(b_sources))
        if source_a == source_b:
            return None
        if any(
            assignments.get((source_a, player)) != node
            for player in coalition_members(block_a, n)
        ) or any(
            assignments.get((source_b, player)) != node
            for player in coalition_members(block_b, n)
        ):
            return None
        relations.append(
            {"source_a": source_a, "source_b": source_b, "sink": node}
        )

    successors: dict[int, int] = {}
    predecessors: dict[int, int] = {}
    relation_by_pair: dict[tuple[int, int], int] = {}
    for relation in relations:
        source_a = relation["source_a"]
        source_b = relation["source_b"]
        if source_a in successors or source_b in predecessors:
            return None
        successors[source_a] = source_b
        predecessors[source_b] = source_a
        relation_by_pair[source_a, source_b] = relation["sink"]
    first = [source for source in sources if source not in predecessors]
    last = [source for source in sources if source not in successors]
    if len(first) != 1 or len(last) != 1:
        return None
    order = [first[0]]
    while order[-1] in successors:
        order.append(successors[order[-1]])
        if len(order) > len(sources):
            return None
    if len(order) != len(sources) or order[-1] != last[0]:
        return None
    if len(relations) != len(sources) - 1:
        return None

    if any(
        assignments.get((order[0], player)) != int(sink_b["node"])
        for player in coalition_members(block_b, n)
    ) or any(
        assignments.get((order[-1], player)) != int(sink_a["node"])
        for player in coalition_members(block_a, n)
    ):
        return None

    expected = set()
    for index, source in enumerate(order):
        if index == 0:
            for player in coalition_members(block_b, n):
                expected.add((source, player, int(sink_b["node"])))
        else:
            previous_sink = relation_by_pair[order[index - 1], source]
            for player in coalition_members(block_b, n):
                expected.add((source, player, previous_sink))
        if index == len(order) - 1:
            for player in coalition_members(block_a, n):
                expected.add((source, player, int(sink_a["node"])))
        else:
            next_sink = relation_by_pair[source, order[index + 1]]
            for player in coalition_members(block_a, n):
                expected.add((source, player, next_sink))
    actual = {
        (source, player, sink)
        for (source, player), sink in assignments.items()
    }
    if actual != expected:
        return None
    return {
        "block_a": block_a,
        "block_b": block_b,
        "source_order": order,
        "left_sink_b": int(sink_b["node"]),
        "right_sink_a": int(sink_a["node"]),
        "intermediate_sinks": [
            relation_by_pair[order[index], order[index + 1]]
            for index in range(len(order) - 1)
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    archive = json.loads(args.archive.read_text())
    report = json.loads(args.terminal_report.read_text())
    n = int(archive["n"])
    grand = (1 << n) - 1
    terminals = report["terminals"]
    sources = [
        int(terminal["node"])
        for terminal in terminals
        if F(terminal["coefficient"]) > 0
    ]
    sinks = [
        {
            "node": int(terminal["node"]),
            "coalition": int(terminal["coalition"]),
            "coefficient": -F(terminal["coefficient"]),
        }
        for terminal in terminals
        if F(terminal["coefficient"]) < 0
    ]
    proper = [sink for sink in sinks if sink["coalition"] != grand]
    if (
        len(proper) != 2
        or proper[0]["coalition"] ^ proper[1]["coalition"] != grand
        or any(sink["coefficient"] != 1 for sink in sinks)
    ):
        raise RuntimeError("dual does not have two complementary unit endpoints")

    path_mass: dict[tuple[int, int, int], F] = defaultdict(F)
    for path in report["terminal_path_decomposition"]:
        path_mass[
            int(path["source"]), int(path["player"]), int(path["sink"])
        ] += F(path["scaled_flow"])
    assignments: dict[tuple[int, int], int] = {}
    for source in sources:
        for player in range(n):
            destinations = [
                sink
                for (candidate, candidate_player, sink), mass in path_mass.items()
                if candidate == source and candidate_player == player and mass
            ]
            if len(destinations) != 1:
                raise RuntimeError("flow is not a unit player-path decomposition")
            sink = destinations[0]
            if path_mass[source, player, sink] != 1:
                raise RuntimeError("a player path has non-unit scaled flow")
            assignments[source, player] = sink

    orientation = None
    for block_a, block_b in (
        (proper[0]["coalition"], proper[1]["coalition"]),
        (proper[1]["coalition"], proper[0]["coalition"]),
    ):
        orientation = recognize_orientation(
            sources, sinks, assignments, n, block_a, block_b
        )
        if orientation is not None:
            break
    if orientation is None:
        raise RuntimeError("terminal flow is not an alternating ladder")

    games = [
        [F(value) for value in game] for game in archive["family"]["games"]
    ]
    grand_values = {game[grand] for game in games}
    if len(grand_values) != 1:
        raise RuntimeError("alternating-ladder theorem requires a common grand")
    common_grand = next(iter(grand_values))
    block_a = orientation["block_a"]
    block_b = orientation["block_b"]
    order = orientation["source_order"]
    left = orientation["left_sink_b"]
    right = orientation["right_sink_a"]
    intermediate = orientation["intermediate_sinks"]

    checks = []
    for index, sink in enumerate(intermediate):
        source_a = order[index]
        source_b = order[index + 1]
        checks.append(
            {
                "source_a": source_a,
                "source_b": source_b,
                "sink": sink,
                "sink_b_equals_source_a_b": (
                    games[sink][block_b] == games[source_a][block_b]
                ),
                "sink_a_equals_source_b_a": (
                    games[sink][block_a] == games[source_b][block_a]
                ),
                "source_b_b_le_source_a_b": (
                    games[source_b][block_b] <= games[source_a][block_b]
                ),
            }
        )
    check_fields = (
        "sink_b_equals_source_a_b",
        "sink_a_equals_source_b_a",
        "source_b_b_le_source_a_b",
    )
    if not all(
        all(bool(check[field]) for field in check_fields) for check in checks
    ):
        raise RuntimeError("an exact protected-path ladder identity failed")
    if games[left][block_a] != games[order[0]][block_a]:
        raise RuntimeError("left endpoint invariance failed")
    if games[right][block_b] != games[order[-1]][block_b]:
        raise RuntimeError("right endpoint invariance failed")

    residual = games[left][block_a] + games[right][block_b] - common_grand
    first_balanced_slack = (
        common_grand
        - games[order[0]][block_a]
        - games[order[0]][block_b]
    )
    monotonic_slack = games[order[0]][block_b] - games[order[-1]][block_b]
    dual_residual = F(report["dual_objective"]) * int(report["scale"])
    if residual != dual_residual or residual > 0:
        raise RuntimeError("ladder objective is not exactly certified nonpositive")
    if -residual != first_balanced_slack + monotonic_slack:
        raise RuntimeError("ladder slack decomposition failed")

    result = {
        "status": "exact_alternating_ladder_audit",
        "source": str(args.archive),
        "terminal_report": str(args.terminal_report),
        "theorem_applies": True,
        "common_grand": str(common_grand),
        "block_a": block_a,
        "block_b": block_b,
        "source_order": order,
        "left_sink_b": left,
        "right_sink_a": right,
        "intermediate_sinks": intermediate,
        "dual_scale": int(report["scale"]),
        "dual_objective_unnormalized": str(residual),
        "safe_slack": str(-residual),
        "safe_slack_decomposition": {
            "first_source_balanced_slack": str(first_balanced_slack),
            "block_b_monotonic_slack": str(monotonic_slack),
        },
        "protected_path_checks": checks,
        "interpretation": (
            "Protected-path invariance makes block B nonincreasing along the "
            "source ladder. The endpoint objective is therefore bounded by "
            "the first source's A/B balancedness inequality."
        ),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
