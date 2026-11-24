"""Smoke-test the GTDB preparation shell using a local MMseqs2 stub."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent


class PrepChimeraBuildTests(unittest.TestCase):
    def test_build_uses_local_inputs_custom_prefix_without_wget(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = root / 'downloaded GTDB R232'
            inputs.mkdir()
            (inputs / 'gtdb_proteins_aa_reps_r232.tar.gz').write_bytes(b'fake archive')
            (inputs / 'bac120_taxonomy_r232.tsv').write_text(
                'GB_GCA_0001\td__Bacteria;p__Nitrospirota;c__Nitrospiria;o__Nitrospirales;'
                'f__Nitrospiraceae;g__Nitrospira;s__Nitrospira sp001\n')
            (inputs / 'ar53_taxonomy_r232.tsv').write_text(
                'RS_GCF_0002\td__Archaea;p__Thermoproteota;c__Nitrososphaeria;'
                'o__Nitrososphaerales;f__Nitrososphaeraceae;g__Nitrososphaera;s__Nitrososphaera sp002\n')
            bindir = root / 'bin'
            bindir.mkdir()
            exe = bindir / 'mmseqs'
            exe.write_text('''#!/usr/bin/env bash
set -euo pipefail
case "$1" in
  tar2db)
    [[ -s "$2" ]]
    printf '1\\tGB_GCA_0001_protein.faa.gz\\n2\\tRS_GCF_0002_protein.faa.gz\\n' > "$3.lookup"
    ;;
  createdb)
    [[ -f "$2.lookup" ]]
    printf 'ok' > "$3.dbtype"
    ;;
  createtaxdb)
    [[ -f "$2.dbtype" ]]
    printf 'ok' > "${2}_mapping"
    printf 'ok' > "${2}_taxonomy"
    ;;
  *) exit 10 ;;
esac
''')
            exe.chmod(0o755)
            wget = bindir / 'wget'
            wget.write_text('#!/bin/sh\nexit 99\n')
            wget.chmod(0o755)
            prefix = root / 'shared db' / 'GTDB_R232'
            result = subprocess.run(['bash', str(ROOT / 'modules/00_prep.sh'), 'chimera',
                                     str(root / 'work'), '1', str(inputs), str(prefix)],
                                    env=dict(os.environ, PATH=str(bindir)+os.pathsep+os.environ['PATH']),
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertEqual((Path(str(prefix)+'.version')).read_text(), 'GTDB release R232 (232.0)\n')
            self.assertTrue((Path(str(prefix)+'_taxonomy')).is_file())
            self.assertTrue((Path(str(prefix)+'_mapping')).is_file())
            lookup = (root / 'work/tardb.lookup').read_text()
            self.assertNotIn('_protein.faa.gz', lookup)
            mapping = (root / 'work/taxonomy/mapping_genomes').read_text().splitlines()
            self.assertEqual(len(mapping), 2)
            self.assertTrue(mapping[0].startswith('GB_GCA_0001\t'))
            self.assertTrue(mapping[1].startswith('RS_GCF_0002\t'))
            self.assertTrue(all((inputs / f).is_file() for f in (
                'gtdb_proteins_aa_reps_r232.tar.gz', 'bac120_taxonomy_r232.tsv', 'ar53_taxonomy_r232.tsv')))


if __name__ == '__main__':
    unittest.main()
