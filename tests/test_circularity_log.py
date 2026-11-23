"""Regression tests for the circularity module's shared run log."""
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / 'modules' / '02_circularity.sh'


class CircularityLogTests(unittest.TestCase):
    def test_run_creates_log_and_replaces_previous_run(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / 'output'
            env = dict(os.environ, OUT=str(output),
                       MAG_DIR=str(Path(temp) / 'nonexistent_mags'),
                       READS=str(Path(temp) / 'nonexistent.fastq'))
            command = ['bash', str(SCRIPT), '--outdir', str(output / '02_circularity')]
            log = output / '00_log' / '02_circularity.log'

            result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=20)
            self.assertNotEqual(result.returncode, 0)
            self.assertTrue(log.is_file())
            self.assertIn('Directory not found:', log.read_text())

            log.write_text('old-run-marker\n')
            result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=20)
            self.assertNotEqual(result.returncode, 0)
            self.assertNotIn('old-run-marker', log.read_text())

    def test_worker_does_not_truncate_parent_log(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / 'output'
            log = output / '00_log' / '02_circularity.log'
            log.parent.mkdir(parents=True)
            log.write_text('parent-run-marker\n')
            env = dict(os.environ, OUT=str(output),
                       MAG_DIR=str(Path(temp) / 'nonexistent_mags'),
                       READS=str(Path(temp) / 'nonexistent.fastq'))
            result = subprocess.run(
                ['bash', str(SCRIPT), '--group-worker', '--outdir', str(output / '02_circularity')],
                env=env, capture_output=True, text=True, timeout=20)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(log.read_text(), 'parent-run-marker\n')


if __name__ == '__main__':
    unittest.main()
