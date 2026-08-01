#!/usr/bin/env python3
"""Find small rational perspective anchors with feasible local dual blocks."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog

from n5_facet_search import load_facets
from n5_mixed_terminal_branch_search import extreme_representations


F = Fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("terminal_report", type=Path)
    parser.add_argument("float_support", type=Path)
    parser.add_argument("--denominators", default="100,1000,10000,100000")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = json.loads(args.terminal_report.read_text())
    support = json.loads(args.float_support.read_text())
    terminals = report["terminals"]
    n = len(terminals[0]["scaled_divergence"])
    grand = (1 << n) - 1
    coalition_count = grand + 1
    facets, _ = load_facets()
    mixed = [
        terminal
        for terminal in terminals
        if terminal.get("coefficient") is None
        and len(set(terminal["scaled_divergence"])) > 2
    ]
    representations = {
        int(terminal["node"]): extreme_representations(
            [F(value) for value in terminal["scaled_divergence"]], n
        )
        for terminal in mixed
    }
    float_simplex: dict[int, float] = {}
    float_links: dict[int, list[float]] = {
        int(terminal["node"]): [0.0] * coalition_count
        for terminal in mixed
    }
    for active in support["equalities"]:
        metadata = active["row"]
        if metadata[0] == "perspective_simplex":
            float_simplex[int(metadata[1])] = float(active["weight"])
        elif metadata[0] == "perspective_link":
            float_links[int(metadata[1])][int(metadata[2])] = float(
                active["weight"]
            )

    facet_count = len(facets)
    empty_variable = facet_count
    grand_variable = facet_count + 1
    choice_nonnegative_variable = facet_count + 2
    local_variable_count = facet_count + 3
    stationarity = np.zeros((coalition_count + 1, local_variable_count))
    for facet_index, facet in enumerate(facets):
        for coalition, coefficient in enumerate(facet):
            stationarity[coalition, facet_index] = -float(coefficient)
    stationarity[0, empty_variable] = 1.0
    stationarity[grand, grand_variable] = 1.0
    stationarity[-1, grand_variable] = -1.0
    stationarity[-1, choice_nonnegative_variable] = -1.0
    bounds = (
        [(None, 0.0)] * facet_count
        + [(None, None), (None, None), (None, 0.0)]
    )

    rows: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    for denominator in (
        int(value) for value in args.denominators.split(",")
    ):
        simplex = {
            node: F(value).limit_denominator(denominator)
            for node, value in float_simplex.items()
        }
        links = {
            node: [
                F(value).limit_denominator(denominator) for value in values
            ]
            for node, values in float_links.items()
        }
        feasible_count = 0
        first_failure = None
        for node, choices in representations.items():
            for choice, (mu, coalitions) in enumerate(choices):
                weights = [F(0)] * coalition_count
                for coalition, weight in coalitions:
                    weights[coalition] = weight
                rhs = np.asarray(
                    [
                        float(-weights[coalition] + links[node][coalition])
                        for coalition in range(coalition_count)
                    ]
                    + [float(-mu - simplex[node])]
                )
                result = linprog(
                    np.zeros(local_variable_count),
                    A_eq=stationarity,
                    b_eq=rhs,
                    bounds=bounds,
                    method="highs-ds",
                    options={
                        "dual_feasibility_tolerance": 1e-10,
                        "primal_feasibility_tolerance": 1e-10,
                    },
                )
                if result.success:
                    feasible_count += 1
                elif first_failure is None:
                    first_failure = {
                        "node": node,
                        "choice": choice,
                        "message": result.message,
                    }
        total = sum(len(choices) for choices in representations.values())
        row = {
            "max_denominator": denominator,
            "feasible_local_blocks": feasible_count,
            "total_local_blocks": total,
            "first_failure": first_failure,
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
        if feasible_count == total:
            best = {
                "max_denominator": denominator,
                "simplex": {
                    str(node): str(value) for node, value in simplex.items()
                },
                "links": {
                    str(node): [str(value) for value in values]
                    for node, values in links.items()
                },
            }
            break

    payload = {
        "status": "perspective_rational_anchor_search",
        "source": str(args.terminal_report),
        "float_support": str(args.float_support),
        "rows": rows,
        "best_complete_anchor": best,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"best_complete_anchor": best}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
