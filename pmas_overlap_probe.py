#!/usr/bin/env python3
"""Probe the overlapping-subgame obstruction of Hokari and Ishida (2026)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def members(mask: int) -> list[int]:
    return [player for player in range(5) if mask >> player & 1]


def matrix(rows: list[dict[int, float]], columns: int):
    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[float] = []
    for row_index, row in enumerate(rows):
        for column, value in row.items():
            if value:
                row_indices.append(row_index)
                column_indices.append(column)
                values.append(value)
    return coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(rows), columns),
    ).tocsr()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    worth = {
        0: 0,
        1: 0,
        2: 0,
        4: 0,
        8: 0,
        16: 0,
        3: 0,
        5: 6,
        9: 0,
        6: 6,
        10: 0,
        12: 6,
        7: 9,
        11: 0,
        13: 9,
        14: 9,
        15: 11,
        20: 0,
        24: 6,
        28: 9,
    }
    hosts = (15, 28)
    coalitions = sorted(
        coalition
        for coalition in worth
        if coalition and any(coalition & host == coalition for host in hosts)
    )
    variables: list[tuple[int, int]] = []
    for coalition in coalitions:
        for player in members(coalition):
            variables.append((coalition, player))
    index = {variable: position for position, variable in enumerate(variables)}
    margin_index = len(variables)
    variable_count = margin_index + 1

    inequalities: list[dict[int, float]] = []
    rhs: list[float] = []
    metadata: list[tuple[Any, ...]] = []
    equalities: list[dict[int, float]] = []
    equality_rhs: list[float] = []
    equality_metadata: list[tuple[Any, ...]] = []
    for coalition in coalitions:
        equalities.append(
            {index[coalition, player]: 1.0 for player in members(coalition)}
        )
        equality_rhs.append(float(worth[coalition]))
        equality_metadata.append(("efficiency", coalition))
        subcoalition = coalition
        while True:
            subcoalition = (subcoalition - 1) & coalition
            if not subcoalition:
                break
            inequalities.append(
                {
                    index[coalition, player]: -1.0
                    for player in members(subcoalition)
                }
            )
            rhs.append(-float(worth[subcoalition]))
            metadata.append(("core", coalition, subcoalition))
    for lower in coalitions:
        for upper in coalitions:
            if lower == upper or lower & upper != lower:
                continue
            if (upper ^ lower).bit_count() != 1:
                continue
            for player in members(lower):
                inequalities.append(
                    {
                        index[lower, player]: 1.0,
                        index[upper, player]: -1.0,
                        margin_index: 1.0,
                    }
                )
                rhs.append(0.0)
                metadata.append(("population_monotonicity", lower, upper, player))
    objective = np.zeros(variable_count)
    objective[margin_index] = -1.0
    solved = linprog(
        objective,
        A_ub=matrix(inequalities, variable_count),
        b_ub=np.asarray(rhs),
        A_eq=matrix(equalities, variable_count),
        b_eq=np.asarray(equality_rhs),
        bounds=[(None, None)] * variable_count,
        method="highs-ds",
    )
    if not solved.success:
        raise RuntimeError(solved.message)
    active_inequalities = [
        {
            "row": list(metadata[row]),
            "dual": float(solved.ineqlin.marginals[row]),
            "slack": float(solved.ineqlin.residual[row]),
        }
        for row in range(len(inequalities))
        if abs(float(solved.ineqlin.marginals[row])) > 1e-9
    ]
    active_equalities = [
        {
            "row": list(equality_metadata[row]),
            "dual": float(solved.eqlin.marginals[row]),
        }
        for row in range(len(equalities))
        if abs(float(solved.eqlin.marginals[row])) > 1e-9
    ]
    result = {
        "status": "hokari_ishida_overlap_pmas_probe",
        "hosts": list(hosts),
        "coalitions": coalitions,
        "maximum_common_margin": float(solved.x[margin_index]),
        "weak_pmas_exists": bool(solved.x[margin_index] >= -1e-9),
        "active_inequalities": active_inequalities,
        "active_equalities": active_equalities,
        "allocations": {
            str(coalition): {
                str(player): float(solved.x[index[coalition, player]])
                for player in members(coalition)
            }
            for coalition in coalitions
        },
    }
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
