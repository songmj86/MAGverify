"""Parallel module loops must finish while stdout is sent to a `tee` log."""
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parent.parent


class BackgroundWaitTests(unittest.TestCase):
    def test_jobs_finish_with_tee_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            script = f'''
set -Eeuo pipefail
source "{ROOT}/modules/lib_qc.sh"
exec > >(tee -a "{tmp}/log.txt") 2>&1
for i in 1 2 3 4 5; do ( sleep 0.3; : > "{tmp}/$i.done" ) & qc_wait_for_slot 2; done
qc_wait_all
ls "{tmp}"/*.done | wc -l
'''
            start = time.monotonic()
            result = subprocess.run(['bash', '-c', script], capture_output=True, text=True, timeout=20)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.strip().splitlines()[-1], '5')
            self.assertGreater(time.monotonic() - start, 0.8)

    def test_modules_do_not_use_bare_wait(self):
        for path in sorted((ROOT / 'modules').glob('*.sh')):
            for number, line in enumerate(path.read_text().splitlines(), 1):
                stripped = line.strip()
                with self.subTest(file=path.name, line=number):
                    self.assertNotIn(stripped, ('wait', 'wait -n', 'wait -n || true'))


if __name__ == '__main__':
    unittest.main()
