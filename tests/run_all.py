"""
tests/run_all.py
================

Run every test in this directory and report pass/fail.
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
TESTS = [
    "test_taxonomy.py",
    "test_step2_pipeline.py",
    "test_step3_pipeline.py",
    "test_step4_pipeline.py",
    "test_pipeline_full.py",
]


def main() -> int:
    passed: list[str] = []
    failed: list[tuple[str, str]] = []

    for t in TESTS:
        path = TESTS_DIR / t
        if not path.exists():
            failed.append((t, "file not found"))
            continue
        print(f"\n{'#' * 78}\n# {t}\n{'#' * 78}")
        proc = subprocess.run([sys.executable, str(path)], cwd=TESTS_DIR.parent)
        if proc.returncode == 0:
            passed.append(t)
        else:
            failed.append((t, f"exited {proc.returncode}"))

    print("\n" + "=" * 78)
    print(f"PASSED: {len(passed)}/{len(TESTS)}")
    for t in passed:
        print(f"  ✓ {t}")
    if failed:
        print(f"\nFAILED: {len(failed)}")
        for t, why in failed:
            print(f"  ✗ {t}  ({why})")
    print("=" * 78)
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
