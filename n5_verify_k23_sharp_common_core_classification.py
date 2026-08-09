#!/usr/bin/env python3
"""Verify exact common-core coverage of every sharp five-player K2,3 orbit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from n5_k23_batch import canonical_topology, exhaustive_orbits


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=Path("results"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    representatives, sink_normalized_count = exhaustive_orbits()
    expected = set(representatives)
    paths = sorted(
        [args.results / "n5_k23_all_pairs_sharp_common_core_exact.json"]
        + list(args.results.glob("n5_k23_orbit_*_common_core_exact.json"))
    )
    records: dict[tuple[int, ...], dict[str, object]] = {}
    for path in paths:
        if not path.exists():
            continue
        payload = json.loads(path.read_text())
        if payload.get("status") != "k23_common_core_exact":
            raise AssertionError(f"{path} is not an exact common-core result")
        if payload.get("common_core_gap_exact") != "0":
            raise AssertionError(f"{path} has a nonzero common-core gap")
        if payload.get("balanced_vertex_count") != 1_291:
            raise AssertionError(f"{path} has the wrong balanced-vertex count")
        if payload.get("exact_failure_count") != 0:
            raise AssertionError(f"{path} has exact reconstruction failures")
        if payload.get("exact_certificate_count") != payload.get("solved"):
            raise AssertionError(f"{path} does not certify every solved case")
        if payload.get("solved") != payload.get("envelope_assignment_count"):
            raise AssertionError(f"{path} did not exhaust its envelope cases")
        representative = canonical_topology(tuple(payload["masks"]))
        if representative in records:
            raise AssertionError(f"duplicate orbit certificate for {representative}")
        records[representative] = {
            "representative": list(representative),
            "source": str(path),
            "varying_sink_coordinates": payload["varying_sink_coordinates"],
            "envelope_assignment_count": payload["envelope_assignment_count"],
            "exact_reconstruction_max_denominator": payload[
                "exact_reconstruction_max_denominator"
            ],
        }
    observed = set(records)
    if observed != expected:
        raise AssertionError(
            f"orbit coverage mismatch: missing={sorted(expected - observed)}, "
            f"extra={sorted(observed - expected)}"
        )
    cache = json.loads((args.results / "n5_balanced_vertices_exact.json").read_text())
    if cache.get("status") != "n5_balanced_vertices_exact":
        raise AssertionError("balanced-vertex cache has the wrong status")
    if cache.get("vertex_count") != 1_291 or len(cache.get("vertices", [])) != 1_291:
        raise AssertionError("balanced-vertex cache has the wrong size")
    total_cases = sum(
        int(record["envelope_assignment_count"]) for record in records.values()
    )
    payload = {
        "status": "n5_k23_all_sharp_common_core_classification_exact",
        "sink_normalized_configuration_count": sink_normalized_count,
        "symmetry_orbit_count": len(expected),
        "balanced_vertex_count": 1_291,
        "exact_envelope_certificate_count": total_cases,
        "exact_failure_count": 0,
        "all_common_core_gaps_exact": "0",
        "orbits": [records[representative] for representative in representatives],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        json.dumps(
            {
                key: payload[key]
                for key in (
                    "status",
                    "sink_normalized_configuration_count",
                    "symmetry_orbit_count",
                    "balanced_vertex_count",
                    "exact_envelope_certificate_count",
                    "exact_failure_count",
                    "all_common_core_gaps_exact",
                )
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
