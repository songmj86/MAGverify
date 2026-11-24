"""Without -t or -m, MAGverify uses 80% of the detected CPUs and memory."""
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
LIB = ROOT / 'modules' / 'lib_config.sh'


def source_lib(tmp, **env):
    script = (f'source "{LIB}"\n'
              'printf "%s %s %s %s\\n" "$MAX_THREADS" "$THREADS" "$MAX_MEMORY_GB" "$MEMORY_GB"\n')
    return subprocess.run(['bash', '-c', script], capture_output=True, text=True,
                          env={'PATH': '/usr/bin:/bin', 'OUT': str(tmp), **env})


class ResourceDefaultTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.tmp = Path(self.temp.name)

    def values(self, **env):
        result = source_lib(self.tmp, **env)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout.split()

    def test_threads_default_to_eighty_percent(self):
        max_threads, threads, _max_mem, _mem = self.values(SLURM_CPUS_PER_TASK='10')
        self.assertEqual((max_threads, threads), ('10', '8'))

    def test_threads_never_fall_below_one(self):
        _max_threads, threads, _max_mem, _mem = self.values(SLURM_CPUS_PER_TASK='1')
        self.assertEqual(threads, '1')

    def test_explicit_threads_may_use_the_whole_allocation(self):
        _max_threads, threads, _max_mem, _mem = self.values(SLURM_CPUS_PER_TASK='10', THREADS='10')
        self.assertEqual(threads, '10')

    def test_threads_above_the_allocation_are_rejected(self):
        result = source_lib(self.tmp, SLURM_CPUS_PER_TASK='10', THREADS='11')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('exceeds the detected allocation', result.stderr)

    def test_memory_defaults_to_eighty_percent(self):
        _max_threads, _threads, max_mem, mem = self.values()
        if int(max_mem) * 80 // 100 < 1:
            self.skipTest('less than 2 GB of memory detected')
        self.assertEqual(int(mem), int(max_mem) * 80 // 100)

    def test_explicit_memory_above_the_machine_is_rejected(self):
        result = source_lib(self.tmp, MEMORY_GB='999999')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('exceeds the detected memory', result.stderr)

    def test_percentage_is_configurable(self):
        _max_threads, threads, _max_mem, _mem = self.values(SLURM_CPUS_PER_TASK='10',
                                                            RESOURCE_DEFAULT_PERCENT='50')
        self.assertEqual(threads, '5')


if __name__ == '__main__':
    unittest.main()
