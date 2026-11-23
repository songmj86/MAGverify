"""Shared Pyrodigal gene prediction and anvi'o external gene calls."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent

STUB = r'''#!/usr/bin/env bash
if [[ "$1" == "--version" ]]; then echo "pyrodigal v3.6.3"; exit 0; fi
while [[ $# -gt 0 ]]; do case "$1" in -i) in=$2; shift 2;; -a) faa=$2; shift 2;; -o) gff=$2; shift 2;; *) shift;; esac; done
echo call >> "$(dirname "$0")/calls.log"
printf '>ctgA_1 # 1 # 99 # 1 # ID=1_1;partial=00\nMKV*\n>ctgA_2 # 151 # 500 # -1 # ID=1_2;partial=10\nMLL\n>ctgB_1 # 1 # 99 # 1 # ID=2_1;partial=00\nMKV*\n' > "$faa"
printf '##gff-version  3\nctgA\tp\tCDS\t1\t99\t.\t+\t0\tID=1_1;partial=00\nctgA\tp\tCDS\t151\t500\t.\t-\t0\tID=1_2;partial=10\nctgB\tp\tCDS\t1\t99\t.\t+\t0\tID=2_1;partial=00\n' > "$gff"
'''


class GeneCallTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.tmp = Path(self.temp.name)
        bindir = self.tmp / 'bin'
        bindir.mkdir()
        stub = bindir / 'pyrodigal'
        stub.write_text(STUB)
        stub.chmod(0o755)
        self.calls = bindir / 'calls.log'
        self.env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH'])
        self.fasta = self.tmp / 'all.fa'
        self.fasta.write_text('>ctgA\n' + 'A' * 500 + '\n>ctgB\n' + 'C' * 200 + '\n')
        self.genes = self.tmp / 'genes'

    def lib_genes(self, *args):
        return subprocess.run(['bash', str(ROOT / 'modules/lib_genes.sh'), *map(str, args)],
                              text=True, capture_output=True, env=self.env)

    def test_prediction_runs_once_and_is_reused(self):
        for _ in range(2):
            result = self.lib_genes(self.fasta, self.genes, 2)
            self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls.read_text().count('call'), 1)
        self.assertIn('reused', result.stdout)
        rows = (self.genes / 'genes.coordinates.tsv').read_text().splitlines()
        self.assertEqual(rows[0], 'ctgA_1\tctgA\t1\t99')
        expected = self.lib_genes('--fingerprint', self.fasta).stdout.strip()
        self.assertEqual((self.genes / '.cache/genes.inputs.sha256').read_text().strip(), expected)
        self.assertEqual(self.lib_genes('--version').stdout.strip(), '3.6.3')

    def test_anvio_external_gene_calls(self):
        self.assertEqual(self.lib_genes(self.fasta, self.genes, 1).returncode, 0)
        name_map = self.tmp / 'name_map.tsv'
        name_map.write_text('c_000000000001\tctgA\n')
        output = self.tmp / 'calls.tsv'
        result = subprocess.run(
            ['python3', str(ROOT / 'scripts/prep_genes.py'), 'anvio',
             str(self.genes / 'genes.gff'), str(self.genes / 'genes.faa'),
             str(name_map), str(output), '--version', '3.6.3'],
            text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [line.split('\t') for line in output.read_text().splitlines()]
        self.assertEqual(rows[0][:6], ['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial'])
        self.assertEqual(rows[1], ['0', 'c_000000000001', '0', '99', 'f', '0', '1', 'pyrodigal', '3.6.3', 'MKV'])
        self.assertEqual(rows[2], ['1', 'c_000000000001', '150', '500', 'r', '1', '1', 'pyrodigal', '3.6.3', 'MLL'])
        self.assertEqual(len(rows), 3)


if __name__ == '__main__':
    unittest.main()
