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
        args = [str(arg) for arg in args]
        env = dict(self.env)
        if '--dry-run' in args:
            args.remove('--dry-run')
            env['MAGVERIFY_DRY_RUN'] = '1'
        return subprocess.run([str(self.cmd), *args],
                              cwd=self.root, text=True, capture_output=True, env=env)

    def test_help_shows_standard_options(self):
        result = self.run_cli('prep-chimera', '-h')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('-d', result.stdout)
        self.assertNotIn('--database-prefix', result.stdout)
        self.assertIn('(database: 00_prep/GTDB_R232)', result.stdout)
        self.assertIn('  -o ', result.stdout)
        self.assertNotIn('--output-dir', result.stdout)
        self.assertIn('  -t ', result.stdout)
        self.assertNotIn('--threads', result.stdout)
        self.assertNotIn('  -a ', result.stdout)
        self.assertNotIn('--dry-run', result.stdout)

    def test_prep_uses_local_inputs_and_default_database_prefix(self):
        result = self.run_cli('prep-chimera', '-d', self.gtdb,
                              '-o', self.output, '-t', '1')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        argv = (self.output / '00_prep/.cache/prep_chimera_work/invocation.txt').read_text().splitlines()
        self.assertEqual(argv, [str(self.output / '00_prep/.cache/prep_chimera_work'), '1', str(self.gtdb),
                                str(self.output / '00_prep/GTDB_R232')])
        self.assertFalse((self.output / '00_prep/.cache/prep_chimera_work/0_prep').exists())

    def test_database_prefix_option_removed(self):
        result = self.run_cli('prep-chimera', '-d', self.gtdb,
                              '--database-prefix', self.output / 'x', '-o', self.output, '-t', '1')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('unrecognized arguments', result.stderr)

    def test_dry_run_no_output_mutation(self):
        result = self.run_cli('prep-chimera', '-d', self.gtdb,
                              '-o', self.output, '-t', '1', '--dry-run')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('00_prep_chimera.sh', result.stdout)
        self.assertIn(str(self.gtdb), result.stdout)
        self.assertFalse(self.output.exists())

    def test_reject_missing_files(self):
        (self.gtdb / FILES[0]).unlink()
        result = self.run_cli('prep-chimera', '-d', self.gtdb,
                              '-o', self.output, '-t', '1')
        self.assertEqual(result.returncode, 2)
        self.assertIn('gtdb_proteins_aa_reps_r232.tar.gz', result.stderr)
        self.assertFalse(self.output.exists())

    def test_reject_missing_gtdb_dir(self):
        result = self.run_cli('prep-chimera', '-o', self.output, '-t', '1')
        self.assertEqual(result.returncode, 2)
        self.assertIn('required: -d', result.stderr)

    def test_reject_bad_threads(self):
        result = self.run_cli('prep-chimera', '-d', self.gtdb,
                              '-o', self.output, '-t', '0')
        self.assertEqual(result.returncode, 2)
        self.assertIn('positive integer', result.stderr)

    def test_prep_appears_separately_in_list(self):
        result = self.run_cli('list')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Preparation:', result.stdout)
        self.assertIn('  prep-chimera', result.stdout)
        self.assertIn('Analysis modules:', result.stdout)
        self.assertIn('  end-to-end ', result.stdout)


if __name__ == '__main__':
    unittest.main()
