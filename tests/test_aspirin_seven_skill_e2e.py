from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
RUNNER = (
    REPOSITORY_ROOT
    / "examples"
    / "aspirin-seven-skill-e2e"
    / "run_case.py"
)


class AspirinSevenSkillE2ETests(unittest.TestCase):
    def test_offline_case_passes_all_handoffs_and_repeatability_checks(self):
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory) / "acceptance"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER),
                    "--output-dir",
                    str(output_dir),
                ],
                cwd=REPOSITORY_ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads(
                (output_dir / "gold_report.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(report["status"], "passed")
            self.assertEqual(len(report["executed_skills"]), 7)
            self.assertEqual(report["run_count"], 2)
            self.assertEqual(report["validators_passed_per_run"], 8)
            self.assertTrue(report["repeatability"]["passed"])
            self.assertTrue(
                all(report["repeatability"]["by_skill"].values())
            )
            self.assertFalse(report["network"]["used"])
            self.assertFalse(any(report["fees"].values()))


if __name__ == "__main__":
    unittest.main()
