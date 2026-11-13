"""CLI tests for GTDB R232 preparation with manually downloaded inputs."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
FILES = ('gtdb_proteins_aa_reps_r232.tar.gz',
         'bac120_taxonomy_r232.tsv', 'ar53_taxonomy_r232.tsv')


class PrepChimeraCLITests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'checkout'
        for path in ('bin', 'scripts', 'modules'):
            (self.root / path).mkdir(parents=True)
        for name in ('bin/MAGverify', 'magverify_cli.py',
                     'scripts/__init__.py', 'scripts/core_magverify_cli.py'):
            shutil.copy2(ROOT / name, self.root / name)
        (self.root / 'bin/MAGverify').chmod(0o755)
        helper = self.root / 'modules/00_prep_chimera.sh'
        helper.write_text('#!/usr/bin/env bash\nset -euo pipefail\n[[ $# == 4 ]]\nmkdir -p "$1"\nprintf \'%s\\n\' "$1" "$2" "$3" "$4" > "$1/invocation.txt"\n')
        self.gtdb = Path(self.temp.name) / 'downloaded GTDB'
        self.gtdb.mkdir()
        for name in FILES:
            (self.gtdb / name).write_text('mock input\n')
        bindir = self.root / 'test_bin'
        bindir.mkdir()
        for exe in ('mmseqs',):
            shim = bindir / exe
            shim.write_text('#!/bin/sh\nexit 0\n')
            shim.chmod(0o755)
        self.env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH'])
        self.cmd = self.root / 'bin/MAGverify'
        self.output = Path(self.temp.name) / 'test output'

    def run_cli(self, *args):
        return subprocess.run([str(self.cmd), *map(str, args)],
                              cwd=self.root, text=True, capture_output=True, env=self.env)

    def test_help_shows_standard_options(self):
        result = self.run_cli('prep-chimera', '--help')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--gtdb-dir', result.stdout)
        self.assertIn('--database-prefix', result.stdout)
        self.assertIn('--output-dir', result.stdout)
        self.assertIn('--threads', result.stdout)
        self.assertNotIn('--contigs-dir', result.stdout)

    def test_prep_uses_local_inputs_and_default_database_prefix(self):
        result = self.run_cli('prep-chimera', '--gtdb-dir', self.gtdb,
                              '--output-dir', self.output, '--threads', '1')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        argv = (self.output / '0_prep/invocation.txt').read_text().splitlines()
        self.assertEqual(argv, [str(self.output / '0_prep'), '1', str(self.gtdb),
                                str(self.output / '0_prep/GTDB_R232')])
        self.assertFalse((self.output / '0_prep/0_prep').exists())

    def test_custom_database_prefix(self):
        external_db = Path(self.temp.name) / 'external database' / 'GTDB_R232'
        result = self.run_cli('prep-chimera', '--gtdb-dir', self.gtdb,
                              '--database-prefix', external_db,
                              '--output-dir', self.output, '--threads', '1')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        argv = (self.output / '0_prep/invocation.txt').read_text().splitlines()
        self.assertEqual(argv[3], str(external_db))

    def test_dry_run_no_output_mutation(self):
        result = self.run_cli('prep-chimera', '--gtdb-dir', self.gtdb,
                              '--output-dir', self.output, '--threads', '1', '--dry-run')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('00_prep_chimera.sh', result.stdout)
        self.assertIn(str(self.gtdb), result.stdout)
        self.assertFalse(self.output.exists())

    def test_reject_missing_files(self):
        (self.gtdb / FILES[0]).unlink()
        result = self.run_cli('prep-chimera', '--gtdb-dir', self.gtdb,
                              '--output-dir', self.output, '--threads', '1')
        self.assertEqual(result.returncode, 2)
        self.assertIn('gtdb_proteins_aa_reps_r232.tar.gz', result.stderr)
        self.assertFalse(self.output.exists())

    def test_reject_missing_gtdb_dir(self):
        result = self.run_cli('prep-chimera', '--output-dir', self.output, '--threads', '1')
        self.assertEqual(result.returncode, 2)
        self.assertIn('required: --gtdb-dir', result.stderr)

    def test_reject_bad_threads(self):
        result = self.run_cli('prep-chimera', '--gtdb-dir', self.gtdb,
                              '--output-dir', self.output, '--threads', '0')
        self.assertEqual(result.returncode, 2)
        self.assertIn('positive integer', result.stderr)

    def test_prep_appears_separately_in_list(self):
        result = self.run_cli('list')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Preparation (run once):', result.stdout)
        self.assertIn('  prep-chimera', result.stdout)
        self.assertIn('Analysis modules (run separately):', result.stdout)


if __name__ == '__main__':
    unittest.main()
