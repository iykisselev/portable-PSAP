import shutil
import subprocess
import unittest
from pathlib import Path


TEST_DIR = Path(__file__).parent


class RScriptTests(unittest.TestCase):
    def test_genotype_models(self) -> None:
        rscript = shutil.which("Rscript")
        if rscript is None:
            self.skipTest("Rscript is unavailable")

        result = subprocess.run(
            [rscript, str(TEST_DIR / "test_genotype_models.R")],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(
            result.returncode,
            0,
            msg=f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )

    def test_acat(self) -> None:
        rscript = shutil.which("Rscript")
        if rscript is None:
            self.skipTest("Rscript is unavailable")
        result = subprocess.run(
            [rscript, str(TEST_DIR / "test_acat.R")],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(
            result.returncode,
            0,
            msg=f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )


if __name__ == "__main__":
    unittest.main()
