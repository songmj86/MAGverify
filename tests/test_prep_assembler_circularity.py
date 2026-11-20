"""prep-circularity with hifiasm-meta GFA, metaMDBG FASTA and metaFlye assembly_info inputs."""
import gzip
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
CLI = ROOT / 'bin/MAGverify'


class PrepAssemblerCircularityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        tmp = Path(self.temp.name)
        self.gfa = tmp / 'Read.p_ctg.gfa'
        self.gfa.write_text('H\tVN:Z:1.0\nS\ts0.ctg000001l\tACGT\tLN:i:4\nS\ts1.ctg000002c\tACGT\tLN:i:4\n')
        self.mdbg = tmp / 'contigs.fasta.gz'
        with gzip.open(self.mdbg, 'wt') as handle:
            handle.write('>ctg11 length=4 coverage=1.0 circular=yes\nACGT\n'
                         '>ctg12 length=4 coverage=1.0 circular=no\nACGT\n')
        self.flye = tmp / 'assembly_info.txt'
        self.flye.write_text('#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n'
                             'contig_1\t4\t10\tY\tN\t1\t*\t1\n'
                             'contig_2\t4\t10\tN\tN\t1\t*\t2\n')
        self.mags = tmp / 'MAGs'
        self.mags.mkdir()
        for mag, contigs in {'Bin1L': ['s1.ctg000002c'], 'Bin2L': ['s0.ctg000001l'],
                             'Bin3L': ['ctg11'], 'Bin4L': ['contig_1', 'contig_2'],
                             'Bin5L': ['contig_1']}.items():
            (self.mags / f'{mag}.fa').write_text(''.join(f'>{c}\nACGT\n' for c in contigs))
        self.out = tmp / 'out'

    def run_cli(self, *sources):
        args = [str(CLI), 'prep-circularity']
        for source in sources:
            args += ['-c', str(source)]
        args += ['-g', str(self.mags), '-o', str(self.out)]
        return subprocess.run(args, cwd=ROOT, text=True, capture_output=True, env=dict(os.environ))

    def test_builds_mag_and_contig_tables(self):
        result = self.run_cli(self.gfa, self.mdbg, self.flye)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        rows = {}
        lines = (self.out / '00_prep/assembler_circularity.MAG.tsv').read_text().splitlines()
        self.assertEqual(lines[0], 'MAG_ID\tAssembler\tAssembler_circularity\tNumber_of_contigs\tContig')
        for line in lines[1:]:
            mag, assembler, call, n, contig = line.split('\t')
            rows[mag] = (call, n)
            if mag == 'Bin3L':
                self.assertEqual(assembler, 'metaMDBG')
        self.assertEqual(rows, {'Bin1L': ('Y', '1'), 'Bin2L': ('N', '1'), 'Bin3L': ('Y', '1'),
                                'Bin4L': ('N', '2'), 'Bin5L': ('Y', '1')})
        contigs = (self.out / '00_prep/assembler_circularity.contig.tsv').read_text()
        self.assertIn('ctg12\tmetaMDBG\tN\tcontigs.fasta.gz', contigs)
        self.assertIn('s1.ctg000002c\thifiasm-meta\tY\tRead.p_ctg.gfa', contigs)
        self.assertIn('contig_2\tmetaFlye\tN\tassembly_info.txt', contigs)

    def test_incomplete_sources_are_rejected(self):
        partial = Path(self.temp.name) / 'headers.txt'
        partial.write_text('S\ts0.ctg000001l\n')
        result = self.run_cli(partial, self.mdbg, self.flye)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('absent from the circularity sources', result.stderr + result.stdout)
        self.assertIn('Bin1L:s1.ctg000002c', result.stderr + result.stdout)


if __name__ == '__main__':
    unittest.main()
