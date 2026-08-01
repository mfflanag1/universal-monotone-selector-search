#!/usr/bin/env python3
"""Classify exact active-component supports in the all-m P/R/T templates."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from fractions import Fraction
from pathlib import Path

from n5_full_m_recurrence_audit import (
    branch_geometry,
    checkpoint_branch_states,
    load,
)


F = Fraction
GRAND = 31


def infer_m(path: Path) -> int:
    match = re.search(r"grid(\d+)", path.name)
    if match is None:
        raise ValueError(f"cannot infer m from {path}")
    return int(match.group(1))


def state_for_game(
    game: tuple[F, ...],
    m: int,
    base: tuple[F, ...],
    endpoints: list[tuple[F, ...]],
    differences: list[list[int]],
) -> tuple[int, tuple[int, ...]]:
    if all(game[coalition] == base[coalition] for coalition in range(1, GRAND)):
        return 0, (0,) * 10
    for branch in (0, 1):
        state = tuple(
            int(
                (game[coalition] - base[coalition])
                * m
                / (endpoints[branch][coalition] - base[coalition])
            )
            for coalition in differences[branch]
        )
        reconstructed = list(base)
        for level, coalition in zip(state, differences[branch], strict=True):
            reconstructed[coalition] += (
                endpoints[branch][coalition] - base[coalition]
            ) * F(level, m)
        if all(
            reconstructed[coalition] == game[coalition]
            for coalition in range(1, GRAND)
        ):
            return branch + 1, state
    raise ValueError("game is not on either intrinsic branch")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m5", type=Path, required=True)
    parser.add_argument("--m6", type=Path, required=True)
    parser.add_argument("--archive", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    checkpoints = load([args.m5, args.m6])
    base, endpoints, differences = branch_geometry()
    templates = []
    for branch in (0, 1):
        states5 = checkpoint_branch_states(
            checkpoints[5][1], 5, branch, base, endpoints, differences
        )
        states6 = checkpoint_branch_states(
            checkpoints[6][1], 6, branch, base, endpoints, differences
        )
        templates.append(
            {
                "P": {state for state in states5 if sum(state) < 24},
                "R": {state for state in states6 if 24 <= sum(state) < 34},
                "T": {state for state in states5 if sum(state) >= 24},
            }
        )

    failures: list[str] = []
    rows = []
    r_support_sets: dict[tuple[int, int, int], set[tuple[int, ...]]] = {}
    for path in sorted(args.archive, key=infer_m):
        m = infer_m(path)
        payload = json.loads(path.read_text())
        games = [tuple(F(value) for value in game) for game in payload["family"]["games"]]
        counts: Counter[str] = Counter()
        branch_counts: dict[str, Counter[str]] = {
            "1": Counter(),
            "2": Counter(),
        }
        for game in games:
            branch, state = state_for_game(
                game, m, base, endpoints, differences
            )
            if branch == 0:
                counts["origin"] += 1
                continue
            template = templates[branch - 1]
            label = None
            base_state = state
            if state in template["P"]:
                label = "P"
            else:
                terminal_base = tuple(value - (m - 5) for value in state)
                if terminal_base in template["T"]:
                    label = "T"
                    base_state = terminal_base
                else:
                    for block in range(m - 5):
                        candidate = tuple(value - block for value in state)
                        if candidate in template["R"]:
                            label = f"R{block}"
                            base_state = candidate
                            r_support_sets.setdefault(
                                (m, branch, block), set()
                            ).add(base_state)
                            break
            if label is None:
                failures.append(
                    f"m={m} branch={branch} unclassified state {state}"
                )
                label = "unclassified"
            counts[label] += 1
            branch_counts[str(branch)][label] += 1
        rows.append(
            {
                "m": m,
                "source": str(path),
                "node_count": len(games),
                "counts": dict(sorted(counts.items())),
                "branch_counts": {
                    branch: dict(sorted(counter.items()))
                    for branch, counter in branch_counts.items()
                },
            }
        )

    block_comparisons = []
    resolutions = sorted({key[0] for key in r_support_sets})
    for branch in (1, 2):
        for left_m, right_m in zip(resolutions, resolutions[1:], strict=False):
            left_last = left_m - 6
            right_last = right_m - 6
            for relation, left_block, right_block in (
                ("fixed_first", 0, 0),
                ("fixed_last", left_last, right_last),
                ("old_last_to_penultimate", left_last, max(0, right_last - 1)),
            ):
                left = r_support_sets.get((left_m, branch, left_block), set())
                right = r_support_sets.get((right_m, branch, right_block), set())
                block_comparisons.append(
                    {
                        "branch": branch,
                        "relation": relation,
                        "left": {"m": left_m, "block": left_block, "count": len(left)},
                        "right": {"m": right_m, "block": right_block, "count": len(right)},
                        "intersection": len(left & right),
                        "left_only": len(left - right),
                        "right_only": len(right - left),
                    }
                )

    result = {
        "status": "pass" if not failures else "fail",
        "rows": rows,
        "repeated_block_support_comparisons": block_comparisons,
        "failures": failures,
        "interpretation": (
            "This classifies finite optimal support archives. Equality of "
            "block supports is evidence for a parametric basis, not itself "
            "an all-m optimality proof."
        ),
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
