"""CLI tests for anvi'o SCG preparation."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent


class PrepScgCLITests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        tmp = Path(self.temp.name)
        self.root = tmp / 'checkout'
        for path in ('bin', 'scripts', 'modules'):
            (self.root / path).mkdir(parents=True)
        for name in ('bin/MAGverify', 'magverify_cli.py',
                     'scripts/__init__.py', 'scripts/core_magverify_cli.py'):
            shutil.copy2(ROOT / name, self.root / name)
        (self.root / 'bin/MAGverify').chmod(0o755)
        (self.root / 'modules/00_prep.sh').write_text(
            '#!/usr/bin/env bash\nset -euo pipefail\n[[ "$1" == scg ]]\nshift\nmkdir -p "$2"\n'
            'printf "%s\\n" "$1" "$3" > "$2/invocation.txt"\n'
            'printf "%s\\n" "$GENE_CALLS_DIR" >> "$2/invocation.txt"\n'
            'printf "contig\\thmm_name\\tgene_callers_id\\tsource\\n" > "$SCG_TSV"\n')
        bindir = self.root / 'test_bin'
        bindir.mkdir()
        for exe in ('pyrodigal', 'anvi-script-reformat-fasta', 'anvi-gen-contigs-database', 'anvi-run-hmms'):
            shim = bindir / exe
            shim.write_text('#!/bin/sh\nexit 0\n')
            shim.chmod(0o755)
        self.env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH'])
        self.assemblies = tmp / 'assemblies'
        self.assemblies.mkdir()
        (self.assemblies / 'asmA.fa').write_text('>a1\nACGT\n')
        (self.assemblies / 'asmB.fasta').write_text('>b1\nGGCC\n')
        self.output = tmp / 'out'

    def run_cli(self, *args):
        return subprocess.run([str(self.root / 'bin/MAGverify'), *map(str, args)],
                              cwd=self.root, text=True, capture_output=True, env=self.env)

    def test_help(self):
        result = self.run_cli('prep-scg', '-h')
        self.assertEqual(result.returncode, 0, result.stderr)
        for option in ('  -a ', '  -t ', '  -o ', '00_prep/scg.contig.tsv'):
            self.assertIn(option, result.stdout)

    def test_builds_all_assemblies_and_scg_table(self):
        result = self.run_cli('prep-scg', '-a', self.assemblies, '-o', self.output, '-t', '1')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        combined = self.output / '00_prep/all_assemblies.fa'
        self.assertIn('>a1', combined.read_text())
        self.assertIn('>b1', combined.read_text())
        invocation = (self.output / '00_prep/.cache/scg/invocation.txt').read_text().split()
        self.assertEqual(invocation, [str(combined), '1', str(self.output / '00_prep')])
        self.assertTrue((self.output / '00_prep/scg.contig.tsv').is_file())
        self.assertFalse((self.output / '00_prep/SCG').exists())

    def test_requires_assemblies(self):
        result = self.run_cli('prep-scg', '-o', self.output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('-a', result.stderr)


if __name__ == '__main__':
    unittest.main()
