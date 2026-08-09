#!/usr/bin/env python3
"""Build independently verified unequal-grand certificates for all K2,3 orbits."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
CLASSIFICATION = RESULTS / "n5_k23_all_sharp_common_core_classification_exact.json"
DEFAULT_SEED = RESULTS / "n5_k23_orbit_00_variable_grand_envelope_exact.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def representatives() -> list[tuple[int, ...]]:
    payload = json.loads(CLASSIFICATION.read_text())
    rows = payload["orbits"]
    if int(payload["symmetry_orbit_count"]) != len(rows):
        raise RuntimeError("classification orbit count mismatch")
    result = [tuple(int(mask) for mask in row["representative"]) for row in rows]
    if len(result) != 20 or any(len(masks) != 6 for masks in result):
        raise RuntimeError("expected twenty six-mask K2,3 representatives")
    return result


def result_path(orbit: int) -> Path:
    return RESULTS / f"n5_k23_orbit_{orbit:02d}_variable_grand_envelope_exact.json"


def verified_metadata(path: Path, masks: tuple[int, ...]) -> dict[str, Any]:
    subprocess.run(
        [sys.executable, "n5_perspective_envelope_certificate.py", str(path)],
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
    )
    payload = json.loads(path.read_text())
    archived_masks = tuple(
        int(pair["player_mask"]) for pair in payload["protected_pairs"]
    )
    if archived_masks != masks:
        raise RuntimeError(
            f"protected-mask mismatch: expected {masks}, got {archived_masks}"
        )
    if payload["exact_maximum"] != "0":
        raise RuntimeError("certificate does not report exact zero maximum")
    try:
        display_path = str(path.relative_to(ROOT))
    except ValueError:
        display_path = str(path)
    return {
        "path": display_path,
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
        "mixed_branch_count": int(payload["mixed_branch_count"]),
        "topology_projection_refinement_count": int(
            payload["topology_projection_refinement_count"]
        ),
    }


def run_logged(command: list[str], log: Path) -> None:
    with log.open("ab") as stream:
        subprocess.run(
            command,
            cwd=ROOT,
            check=True,
            stdout=stream,
            stderr=subprocess.STDOUT,
        )


def build_orbit(
    orbit: int,
    masks: tuple[int, ...],
    seed: Path,
    force: bool,
    reuse_converged_temporaries: bool,
) -> dict[str, Any]:
    destination = result_path(orbit)
    started = time.monotonic()
    if destination.exists() and not force:
        metadata = verified_metadata(destination, masks)
        return {
            "orbit": orbit,
            "masks": list(masks),
            "disposition": "verified_existing",
            "elapsed_seconds": time.monotonic() - started,
            **metadata,
        }

    temp_root = None
    if reuse_converged_temporaries:
        candidates = sorted(
            Path("/private/tmp").glob(f"n5-k23-orbit-{orbit:02d}.*"),
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
        for candidate in candidates:
            candidate_cegar = candidate / "envelope_cegar.json"
            candidate_report = candidate / "terminal_report.json"
            if not candidate_cegar.exists() or not candidate_report.exists():
                continue
            try:
                candidate_payload = json.loads(candidate_cegar.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            if candidate_payload.get("status") == "n5_k23_envelope_cegar_converged":
                temp_root = candidate
                break
    reused_cegar = temp_root is not None
    if temp_root is None:
        temp_root = Path(
            tempfile.mkdtemp(
                prefix=f"n5-k23-orbit-{orbit:02d}.", dir="/private/tmp"
            )
        )
    report = temp_root / "terminal_report.json"
    cegar = temp_root / "envelope_cegar.json"
    exact = temp_root / "envelope_exact.json"
    log = temp_root / "pipeline.log"
    try:
        if not reused_cegar:
            run_logged(
                [
                    sys.executable,
                    "n5_topology_terminal_report.py",
                    "--shape",
                    "k23",
                    "--masks",
                    *(str(mask) for mask in masks),
                    "--output",
                    str(report),
                ],
                log,
            )
            run_logged(
                [
                    sys.executable,
                    "n5_k23_envelope_cegar.py",
                    str(report),
                    "--seed-reference",
                    str(seed),
                    "--output",
                    str(cegar),
                    "--iterations",
                    "40",
                    "--cuts-per-node",
                    "20",
                    "--tolerance",
                    "1e-8",
                ],
                log,
            )
        run_logged(
            [
                sys.executable,
                "n5_perspective_envelope_certificate.py",
                str(exact),
                "--generate-from-report",
                str(report),
                "--link-support",
                str(cegar),
                "--max-denominator",
                "1000000",
                "--zero-threshold",
                "0",
                "--topology-projection-epsilon",
                "1/100",
                "--max-preseeded-memberships-per-node",
                "2",
            ],
            log,
        )
        verified_metadata(exact, masks)
        staging = destination.with_suffix(".json.pending")
        if staging.exists():
            staging.unlink()
        shutil.copyfile(exact, staging)
        os.replace(staging, destination)
        metadata = verified_metadata(destination, masks)
        return {
            "orbit": orbit,
            "masks": list(masks),
            "disposition": "generated_and_verified",
            "reused_converged_cegar": reused_cegar,
            "elapsed_seconds": time.monotonic() - started,
            "temporary_log": str(log),
            **metadata,
        }
    except Exception:
        print(
            json.dumps(
                {
                    "orbit": orbit,
                    "status": "failed",
                    "temporary_directory": str(temp_root),
                    "log": str(log),
                }
            ),
            flush=True,
        )
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--orbit", type=int, action="append")
    parser.add_argument("--seed", type=Path, default=DEFAULT_SEED)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--reuse-converged-temporaries", action="store_true")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=RESULTS / "n5_k23_all_sharp_variable_grand_exact_manifest.json",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    masks_by_orbit = representatives()
    selected = (
        sorted(set(args.orbit))
        if args.orbit is not None
        else list(range(len(masks_by_orbit)))
    )
    if any(not 0 <= orbit < len(masks_by_orbit) for orbit in selected):
        parser.error("orbit index is out of range")
    records: list[dict[str, Any]] = []
    failures = []
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=args.workers
    ) as executor:
        futures = {
            executor.submit(
                build_orbit,
                orbit,
                masks_by_orbit[orbit],
                args.seed,
                args.force,
                args.reuse_converged_temporaries,
            ): orbit
            for orbit in selected
        }
        for future in concurrent.futures.as_completed(futures):
            orbit = futures[future]
            try:
                record = future.result()
            except Exception as error:
                failures.append({"orbit": orbit, "error": repr(error)})
                continue
            records.append(record)
            print(
                json.dumps(
                    {
                        "orbit": orbit,
                        "status": record["disposition"],
                        "elapsed_seconds": round(record["elapsed_seconds"], 3),
                        "sha256": record["sha256"],
                    }
                ),
                flush=True,
            )
    records.sort(key=lambda record: int(record["orbit"]))
    manifest = {
        "status": (
            "n5_k23_all_sharp_variable_grand_exact_certificates"
            if not failures and len(records) == 20
            else "n5_k23_variable_grand_exact_certificate_batch_partial"
        ),
        "classification": str(CLASSIFICATION.relative_to(ROOT)),
        "symmetry_orbit_count": len(masks_by_orbit),
        "selected_orbit_count": len(selected),
        "verified_certificate_count": len(records),
        "failure_count": len(failures),
        "exact_failure_count": 0 if not failures else None,
        "all_terminal_objective_maxima_exact": "0" if not failures else None,
        "orbits": records,
        "failures": failures,
    }
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({key: manifest[key] for key in (
        "status",
        "verified_certificate_count",
        "failure_count",
    )}, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
