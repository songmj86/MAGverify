"""00_prep files written by earlier versions move to the flat layout."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts import core_magverify_cli as cli  # noqa: E402


class PrepLayoutTests(unittest.TestCase):
    def test_earlier_layout_is_flattened(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prep = root / '00_prep'
            old = {
                'genes/genes.faa': 'faa', 'genes/genes.gff': 'gff',
                'genes/genes.coordinates.tsv': 'coord', 'genes/.cache/genes.inputs.sha256': 'fp',
                'genes/.cache/genes.lock': '',
                'SCG/SCG_from_anvio.original_contig.tsv': 'scg', 'SCG/contigs.db': 'db',
                'SCG/00_SCG.log': 'log',
                'gtdb_taxonomy.MAG.tsv': 'mag',
                'gene_taxonomy/GTDB_R232/gene_taxonomy.mmseqs2.GTDB_R232.all_contigs.tsv': 'all',
                'gene_taxonomy/GTDB_R232/gene_taxonomy.mmseqs2.GTDB_R232.MAG_contigs.tsv': 'magc',
                'all_assemblies.fasta': '>c\nA\n',
            }
            for name, text in old.items():
                (prep / name).parent.mkdir(parents=True, exist_ok=True)
                (prep / name).write_text(text)
            cli.migrate_prep_layout(root)
            expected = {
                'genes.faa': 'faa', 'genes.gff': 'gff', 'genes.coordinates.tsv': 'coord',
                '.cache/genes.inputs.sha256': 'fp',
                'scg.contig.tsv': 'scg', '.cache/scg/contigs.db': 'db', '.cache/scg/prep_scg.log': 'log',
                'gtdbtk_taxonomy.MAG.tsv': 'mag', 'taxonomy.gene.tsv': 'all',
                '.cache/gene_taxonomy/gene_taxonomy.mmseqs2.GTDB_R232.MAG_contigs.tsv': 'magc',
            }
            for name, text in expected.items():
                self.assertEqual((prep / name).read_text(), text, name)
            visible = sorted(p.name for p in prep.iterdir() if not p.name.startswith('.'))
            self.assertEqual(visible, ['all_assemblies.fa', 'genes.coordinates.tsv', 'genes.faa',
                                       'genes.gff', 'gtdbtk_taxonomy.MAG.tsv', 'scg.contig.tsv',
                                       'taxonomy.gene.tsv'])
            cli.migrate_prep_layout(root)
            for name, text in expected.items():
                self.assertEqual((prep / name).read_text(), text, name)


if __name__ == '__main__':
    unittest.main()
