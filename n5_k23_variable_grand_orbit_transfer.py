#!/usr/bin/env python3
"""Transfer an exact sharp K2,3 envelope library across all symmetry orbits."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csc_matrix, hstack, lil_matrix, vstack

from n5_k23_batch import exhaustive_orbits, permute_mask
from n5_middle_split_identity_certificate import (
    exact_generators,
    find_certificate,
)
from n5_mixed_terminal_search import complete_bipartite_rectangular_topology
from n5_perspective_envelope_certificate import (
    embedded_topology_generators,
    reconstruct,
    verify as verify_envelope_archive,
)


F = Fraction
N = 5
COALITION_COUNT = 1 << N


def permute_vector(
    vector: tuple[F, ...], permutation: tuple[int, ...]
) -> tuple[F, ...]:
    result = [F(0)] * N
    for player, value in enumerate(vector):
        result[permutation[player]] = value
    return tuple(result)


def permute_coefficients(
    coefficients: tuple[F, ...], permutation: tuple[int, ...]
) -> tuple[F, ...]:
    result = [F(0)] * COALITION_COUNT
    for coalition, value in enumerate(coefficients):
        result[permute_mask(coalition, permutation)] = value
    return tuple(result)


def mapping_permutations(
    source: tuple[F, ...], target: tuple[F, ...]
) -> list[tuple[int, ...]]:
    return [
        permutation
        for permutation in itertools.permutations(range(N))
        if permute_vector(source, permutation) == target
    ]


def topology_data(
    masks: tuple[int, ...],
) -> tuple[dict[int, tuple[F, ...]], list[dict[str, Any]]]:
    topology = complete_bipartite_rectangular_topology(
        masks, (1,) * 6, 2, 3, N
    )
    divergences = {
        node: tuple(F(value) for value in divergence)
        for node, divergence in enumerate(topology.divergence)
    }
    protected_pairs = [
        {
            "source": int(source),
            "sink": int(sink),
            "players": [
                player for player in range(N) if mask >> player & 1
            ],
            "player_mask": int(mask),
        }
        for source, sink, mask in topology.arcs
    ]
    return divergences, protected_pairs


def load_reference(reference_path: Path) -> dict[str, Any]:
    raw = reference_path.read_bytes()
    reference = json.loads(raw)
    if reference["status"] != "n5_terminal_perspective_envelope_exact_certificate":
        raise RuntimeError("reference is not an exact envelope certificate")
    reference_nodes = [int(node) for node in reference["node_ids"]]
    if reference_nodes != list(range(5)):
        raise RuntimeError("reference must use nodes 0 through 4")
    return {
        "path": reference_path,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "payload": reference,
        "nodes": reference_nodes,
        "protected_pairs": reference["protected_pairs"],
        "divergences": {
            int(node): tuple(F(value) for value in values)
            for node, values in reference["terminal_divergences"].items()
        },
        "envelopes": {
            int(node): tuple(F(value) for value in values)
            for node, values in reference["envelopes"].items()
        },
    }


def find_library_transfer(
    reference_index: int,
    reference: dict[str, Any],
    target_divergences: dict[int, tuple[F, ...]],
    topology_matrix: csc_matrix,
    topology_generators: list[tuple[F, ...]],
    topology_metadata: list[tuple[Any, ...]],
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    options: list[dict[str, Any]] = []
    option_indices_by_target: dict[int, list[int]] = {
        node: [] for node in range(5)
    }
    option_indices_by_donor: dict[int, list[int]] = {
        node: [] for node in range(5)
    }
    for target in range(5):
        donors = (0, 1) if target < 2 else (2, 3, 4)
        seen = set()
        for donor in donors:
            for permutation in mapping_permutations(
                reference["divergences"][donor],
                target_divergences[target],
            ):
                envelope = permute_coefficients(
                    reference["envelopes"][donor], permutation
                )
                key = (donor, envelope)
                if key in seen:
                    continue
                seen.add(key)
                option_index = len(options)
                options.append(
                    {
                        "target": target,
                        "donor": donor,
                        "permutation": permutation,
                        "envelope": envelope,
                    }
                )
                option_indices_by_target[target].append(option_index)
                option_indices_by_donor[donor].append(option_index)
    candidate_space_size = 1
    for target in range(5):
        count = len(option_indices_by_target[target])
        if count == 0:
            return None, {
                "candidate_space_size": 0,
                "option_count": len(options),
                "milp_attempt_count": 0,
            }
        candidate_space_size *= count

    generator_count = topology_matrix.shape[1]
    option_count = len(options)
    variable_count = generator_count + option_count
    option_columns = np.zeros((5 * COALITION_COUNT, option_count))
    for option_index, option in enumerate(options):
        offset = int(option["target"]) * COALITION_COUNT
        option_columns[offset : offset + COALITION_COUNT, option_index] = [
            float(value) for value in option["envelope"]
        ]
    coalition_rows = hstack(
        [topology_matrix, csc_matrix(option_columns)], format="csc"
    )
    assignment_rows = lil_matrix((10, variable_count))
    for target in range(5):
        for option_index in option_indices_by_target[target]:
            assignment_rows[target, generator_count + option_index] = 1.0
    for donor in range(5):
        for option_index in option_indices_by_donor[donor]:
            assignment_rows[5 + donor, generator_count + option_index] = 1.0
    base_matrix = vstack(
        [coalition_rows, assignment_rows.tocsc()], format="csc"
    )
    base_lower = np.concatenate((np.zeros(5 * COALITION_COUNT), np.ones(10)))
    base_upper = base_lower.copy()
    objective = np.concatenate(
        (
            np.asarray(
                [
                    1.0 + (index % 97) * 1e-7
                    for index in range(generator_count)
                ]
            ),
            np.asarray(
                [1e-9 * (index + 1) for index in range(option_count)]
            ),
        )
    )
    lower = np.zeros(variable_count)
    upper = np.concatenate(
        (np.full(generator_count, np.inf), np.ones(option_count))
    )
    integrality = np.concatenate(
        (np.zeros(generator_count), np.ones(option_count))
    )
    no_good_rows: list[csc_matrix] = []
    for attempt in range(20):
        constraint_matrix = base_matrix
        constraint_lower = base_lower
        constraint_upper = base_upper
        if no_good_rows:
            constraint_matrix = vstack(
                [base_matrix, *no_good_rows], format="csc"
            )
            constraint_lower = np.concatenate(
                (base_lower, np.full(len(no_good_rows), -np.inf))
            )
            constraint_upper = np.concatenate(
                (base_upper, np.full(len(no_good_rows), 4.0))
            )
        solved = milp(
            objective,
            integrality=integrality,
            bounds=Bounds(lower, upper),
            constraints=LinearConstraint(
                constraint_matrix, constraint_lower, constraint_upper
            ),
            options={
                "disp": False,
                "presolve": True,
                "time_limit": 300.0,
                "mip_rel_gap": 0.0,
                "primal_feasibility_tolerance": 1e-10,
                "dual_feasibility_tolerance": 1e-10,
            },
        )
        if not solved.success:
            return None, {
                "candidate_space_size": candidate_space_size,
                "option_count": option_count,
                "milp_attempt_count": attempt + 1,
                "milp_status": int(solved.status),
                "milp_message": solved.message,
            }
        selected = [
            index
            for index, value in enumerate(solved.x[generator_count:])
            if value > 0.5
        ]
        if len(selected) != 5:
            raise RuntimeError("transfer MILP did not select five envelopes")
        selected_options = sorted(
            (options[index] for index in selected),
            key=lambda option: int(option["target"]),
        )
        target = tuple(
            -value
            for option in selected_options
            for value in option["envelope"]
        )
        try:
            support = find_certificate(
                target, topology_generators, allow_z3=False
            )
        except RuntimeError:
            cut = lil_matrix((1, variable_count))
            for option_index in selected:
                cut[0, generator_count + option_index] = 1.0
            no_good_rows.append(cut.tocsc())
            continue
        return {
            "source_reference_index": reference_index,
            "donor_nodes": [
                int(option["donor"]) for option in selected_options
            ],
            "player_permutations": [
                list(option["permutation"]) for option in selected_options
            ],
            "candidate_space_size": candidate_space_size,
            "option_count": option_count,
            "milp_attempt_count": attempt + 1,
            "global_topology_certificate": [
                {
                    "generator": list(topology_metadata[index]),
                    "weight": str(weight),
                }
                for index, weight in support
            ],
        }, {
            "candidate_space_size": candidate_space_size,
            "option_count": option_count,
            "milp_attempt_count": attempt + 1,
        }
    return None, {
        "candidate_space_size": candidate_space_size,
        "option_count": option_count,
        "milp_attempt_count": 20,
        "milp_status": "exact_reconstruction_failures",
    }


def generate(reference_paths: list[Path], output: Path) -> None:
    if not reference_paths:
        raise RuntimeError("at least one exact envelope reference is required")
    references = [load_reference(path) for path in reference_paths]
    reference_nodes = list(range(5))
    representatives, sink_normalized_count = exhaustive_orbits()
    records = []
    for orbit, masks in enumerate(representatives):
        target_divergences, protected_pairs = topology_data(masks)
        topology_generators, topology_metadata = embedded_topology_generators(
            reference_nodes, N, protected_pairs
        )
        topology_matrix = csc_matrix(
            np.asarray(
                [
                    [
                        float(generator[row])
                        for generator in topology_generators
                    ]
                    for row in range(len(reference_nodes) * COALITION_COUNT)
                ]
            )
        )
        found = None
        reference_attempts = []
        for reference_index, reference in enumerate(references):
            if (
                reference["divergences"] == target_divergences
                and reference["protected_pairs"] == protected_pairs
            ):
                found = {
                    "masks": list(masks),
                    "source_reference_index": reference_index,
                    "donor_nodes": list(reference_nodes),
                    "player_permutations": [list(range(N)) for _ in reference_nodes],
                    "candidate_space_size": 1,
                    "option_count": 5,
                    "milp_attempt_count": 0,
                    "global_topology_certificate": reference["payload"][
                        "global_topology_certificate"
                    ],
                }
                reference_attempts.append(
                    {
                        "candidate_space_size": 1,
                        "option_count": 5,
                        "milp_attempt_count": 0,
                        "direct_reference_reuse": True,
                    }
                )
                break
            found, attempt_record = find_library_transfer(
                reference_index,
                reference,
                target_divergences,
                topology_matrix,
                topology_generators,
                topology_metadata,
            )
            reference_attempts.append(attempt_record)
            if found is not None:
                break
        if found is None:
            raise RuntimeError(
                f"no exact envelope-library transfer found for orbit {orbit}: "
                f"{masks}; attempts {reference_attempts}"
            )
        found["masks"] = list(masks)
        found["reference_attempts"] = reference_attempts
        records.append(found)
        print(
            json.dumps(
                {
                    "orbit": orbit,
                    "reference": found["source_reference_index"],
                    "candidate_space": found["candidate_space_size"],
                    "milp_attempts": found["milp_attempt_count"],
                    "global_support": len(
                        found["global_topology_certificate"]
                    ),
                }
            ),
            flush=True,
        )
    payload = {
        "status": "n5_k23_all_sharp_variable_grand_envelope_library_transfer_exact",
        "n": N,
        "source_references": [
            {
                "path": str(reference["path"]),
                "sha256": reference["sha256"],
                "mixed_branch_count": reference["payload"]["mixed_branch_count"],
            }
            for reference in references
        ],
        "sink_normalized_configuration_count": sink_normalized_count,
        "symmetry_orbit_count": len(representatives),
        "exact_failure_count": 0,
        "all_terminal_objective_maxima_exact": "0",
        "orbits": records,
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")


def verify(archive: Path) -> None:
    payload = json.loads(archive.read_text())
    if payload["status"] != "n5_k23_all_sharp_variable_grand_envelope_library_transfer_exact":
        raise RuntimeError("unexpected transfer archive status")
    references = []
    for archived_reference in payload["source_references"]:
        reference_path = Path(archived_reference["path"])
        reference = load_reference(reference_path)
        if reference["sha256"] != archived_reference["sha256"]:
            raise RuntimeError("reference archive hash mismatch")
        if int(reference["payload"]["mixed_branch_count"]) != int(
            archived_reference["mixed_branch_count"]
        ):
            raise RuntimeError("reference mixed-branch count mismatch")
        verify_envelope_archive(reference_path)
        references.append(reference)
    representatives, sink_normalized_count = exhaustive_orbits()
    if int(payload["sink_normalized_configuration_count"]) != sink_normalized_count:
        raise RuntimeError("sink-normalized count mismatch")
    if int(payload["symmetry_orbit_count"]) != len(representatives):
        raise RuntimeError("orbit count mismatch")
    records = payload["orbits"]
    if len(records) != len(representatives):
        raise RuntimeError("transfer archive omits an orbit")

    local_generators, _ = exact_generators(N)
    local_generator_set = set(local_generators)
    for permutation in itertools.permutations(range(N)):
        if {
            permute_coefficients(generator, permutation)
            for generator in local_generators
        } != local_generator_set:
            raise RuntimeError("exact monotone game cone is not permutation invariant")

    for orbit, (masks, record) in enumerate(
        zip(representatives, records, strict=True)
    ):
        if tuple(record["masks"]) != masks:
            raise RuntimeError(f"orbit {orbit} mask mismatch")
        target_divergences, protected_pairs = topology_data(masks)
        reference_index = int(record["source_reference_index"])
        if not 0 <= reference_index < len(references):
            raise RuntimeError("invalid source reference index")
        reference = references[reference_index]
        donor_nodes = tuple(int(node) for node in record["donor_nodes"])
        mappings = tuple(
            tuple(int(player) for player in permutation)
            for permutation in record["player_permutations"]
        )
        if len(donor_nodes) != 5 or len(mappings) != 5:
            raise RuntimeError("invalid transfer dimensions")
        envelopes = {}
        for target_node, (donor, permutation) in enumerate(
            zip(donor_nodes, mappings, strict=True)
        ):
            if sorted(permutation) != list(range(N)):
                raise RuntimeError("invalid player permutation")
            if permute_vector(
                reference["divergences"][donor], permutation
            ) != target_divergences[target_node]:
                raise RuntimeError("divergence transfer mismatch")
            envelopes[target_node] = permute_coefficients(
                reference["envelopes"][donor], permutation
            )
        topology_generators, topology_metadata = embedded_topology_generators(
            list(range(5)), N, protected_pairs
        )
        topology_lookup = {
            tuple(metadata): index
            for index, metadata in enumerate(topology_metadata)
        }
        target = tuple(
            -value
            for node in range(5)
            for value in envelopes[node]
        )
        if reconstruct(
            record["global_topology_certificate"],
            topology_lookup,
            topology_generators,
        ) != target:
            raise RuntimeError(f"orbit {orbit} global identity failed")
    if int(payload["exact_failure_count"]) != 0:
        raise RuntimeError("archive reports an exact failure")
    if payload["all_terminal_objective_maxima_exact"] != "0":
        raise RuntimeError("archive does not report exact zero maxima")
    print(
        "PASS: all sharp five-player K2,3 unequal-grand transfers, "
        f"{len(representatives)} symmetry orbits"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("--generate-from", type=Path, action="append")
    args = parser.parse_args()
    if args.generate_from is not None:
        generate(args.generate_from, args.archive)
    verify(args.archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
