"""circularity_finalize.py contig writes one row per assembly contig."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parent.parent / 'scripts' / 'circularity_finalize.py'

MAG_TABLE = (
    'ID\tAssembler_reported_circularity\tCircularity_evaluable\tCircularity_eligibility\t'
    'Read_supported_circular_topology\tHigh_confidence_circularity\tAdditional_circular_candidate\t'
    'Circular\tNumber_of_contigs\tContig_ID\tContig_length_bp\tHIFI_END_START_SUPPORTED\n'
    'BinA\tY\tY\tY_single_contig\tY\tY\tN\tY\t1\tc1\t4\tY\n'
    'BinB\tN\tN\tN_not_single_contig\tNA\tN\tN\tN\t2\tNA\tNA\tNA\n'
)


class CircularityContigTests(unittest.TestCase):
    def test_every_contig_gets_a_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'asm').mkdir()
            (tmp / 'mag').mkdir()
            (tmp / 'asm/flye.fa').write_text('>c1\nAAAA\n>c2\nCC\n>c3\nGGG\n')
            (tmp / 'asm/mdbg.fasta').write_text('>d1\nAAAAAA\n')
            (tmp / 'mag/BinA.fa').write_text('>c1\nAAAA\n')
            (tmp / 'mag/BinB.fa').write_text('>c2\nCC\n>c3\nGGG\n')
            (tmp / 'mag.tsv').write_text(MAG_TABLE)
            (tmp / 'assembler.tsv').write_text(
                'Contig\tAssembler\tAssembler_circularity\tFile\nc2\tmetaFlye\tN\tx\nd1\tmetaMDBG\tY\tx\n')
            output = tmp / 'circularity.contig.tsv'
            subprocess.run([sys.executable, str(SCRIPT), 'contig', '--mag-table', str(tmp / 'mag.tsv'),
                            '--mag-dir', str(tmp / 'mag'), '--assembly-dir', str(tmp / 'asm'),
                            '--assembly-fastas', 'flye.fa mdbg.fasta',
                            '--assembler-circularity', str(tmp / 'assembler.tsv'),
                            '--output', str(output)], check=True, capture_output=True)
            with output.open() as handle:
                rows = {row['Contig']: row for row in csv.DictReader(handle, delimiter='\t')}
            self.assertEqual(sorted(rows), ['c1', 'c2', 'c3', 'd1'])
            self.assertEqual((rows['c1']['MAG_ID'], rows['c1']['Circular'], rows['c1']['HIFI_END_START_SUPPORTED'],
                              rows['c1']['Assembler_reported_circularity']), ('BinA', 'Y', 'Y', 'Y'))
            self.assertEqual((rows['c2']['Circularity_eligibility'], rows['c2']['Circular'],
                              rows['c2']['Assembler_reported_circularity']), ('N_not_single_contig', 'NA', 'N'))
            self.assertEqual((rows['d1']['MAG_ID'], rows['d1']['Circularity_eligibility'],
                              rows['d1']['Assembler_reported_circularity'], rows['d1']['Contig_length_bp']),
                             ('NA', 'N_unbinned', 'Y', '6'))


if __name__ == '__main__':
    unittest.main()
