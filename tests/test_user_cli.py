"""Command-line tests using small FASTAs and stub analysis stages."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
STAGES = ('map', 'clipping', 'chimera', 'circularity',
          'false-closure', 'variants', 'phantom', 'repeats', 'false-duplication',
          'summary')


class CommandTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.home = Path(self.tmp.name)
        self.repo = self.home / 'repo'
        (self.repo / 'modules').mkdir(parents=True)
        (self.repo / 'scripts').mkdir(parents=True)
        (self.repo / 'bin').mkdir(parents=True)
        for src in ('bin/MAGverify', 'magverify_cli.py', 'modules/__init__.py',
                    'scripts/__init__.py', 'scripts/core_magverify_cli.py'):
            dst = self.repo / src
            shutil.copy2(ROOT / src, dst)
        self.data = self.home / 'different_locations'
        self.data.mkdir()
        self.assembly = self.data / 'full.fasta'
        self.assembly.write_text('>ctg1\nACGTACGT\n')
        self.source_dir = self.home / 'source_assemblies'
        self.source_dir.mkdir()
        self.source = self.source_dir / 'assembler_one.fa'
        self.source.write_text('>ctg1\nACGTACGT\n')
        self.mag = self.home / 'Bin13L.fna'
        self.mag.write_text('>ctg1\nACGTACGT\n')
        self.reads = self.home / 'reads.fq'
        self.reads.write_text('@r1\nACGT\n+\nIIII\n')
        self.output = self.home / 'output'
        self.map_table = self.home / 'MAG.tsv'
        self.map_table.write_text('MAG_ID\tContig\nBin13L\tctg1\n')
        self.generic_tsv = self.home / 'some.tsv'
        self.generic_tsv.write_text('ID\tcall\nctg1\tY\n')
        self.prefix = self.home / 'db'
        (self.home / 'db.dbtype').write_text('db')
        bam = self.home / 'sample.bam'
        bam.write_text('dummy')
        (self.home / 'sample.bam.inputs.sha256').write_text('provenance')
        self.manifest = self.home / 'mapping.tsv'
        with self.manifest.open('w', encoding='utf-8') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['Assembler','Assembly_FASTA','HiFi_all_BAM','HiFi_primary_BAM',
                             'HiFi_secondary_BAM','HiFi_supplementary_BAM',
                             'HiFi_primary_supplementary_BAM','Output_directory'])
            writer.writerow(['assembler_one',str(self.source),str(bam),str(bam),str(bam),
                             str(bam),str(bam),str(self.home)])
        bin_dir = self.home / 'bin'
        bin_dir.mkdir()
        for program in ('samtools', 'minimap2', 'anvi-script-find-misassemblies', 'sha256sum'):
            stub = bin_dir / program
            stub.write_text('#!/bin/sh\nexit 0\n')
            stub.chmod(0o755)
        self.env = dict(os.environ, PATH=str(bin_dir) + os.pathsep + os.environ['PATH'])
        (self.repo / 'modules/01_mapping.sh').write_text('''#!/bin/bash
set -euo pipefail
[[ -n "$ASSEMBLER_FASTAS" ]] && [[ "$QC_SKIP_LOCAL_CONFIG" == 1 ]]
mkdir -p "$2"
bam="$2/example.bam"
printf 'BAM' > "$bam"
printf 'hash' > "$bam.inputs.sha256"
printf 'Assembler\\tAssembly_FASTA\\tHiFi_all_BAM\\tHiFi_primary_BAM\\tHiFi_secondary_BAM\\tHiFi_supplementary_BAM\\tHiFi_primary_supplementary_BAM\\tOutput_directory\\n' > "$2/assembler_mapping_outputs.tsv"
for name in $ASSEMBLER_FASTAS; do
  printf 'assembler\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "$ASSEMBLER_ASSEMBLY_DIR/$name" "$bam" "$bam" "$bam" "$bam" "$bam" "$2" >> "$2/assembler_mapping_outputs.tsv"
done
''')
        (self.repo / 'modules/03_clipping.sh').write_text('''#!/bin/bash
set -euo pipefail
[[ -s "$ASSEMBLER_MAPPING_MANIFEST" ]]
[[ -f "$MAG_DIR/Bin13L.fna" ]]
[[ "$ASSEMBLY_FA" == "$1" ]]
[[ "$OUT" == "$(dirname "$2")" ]]
mkdir -p "$2"
printf '%s\\n%s\\n%s\\n' "${READS:-}" "$MAG_DIR" "${ASSEMBLER_FASTAS:-}" > "$2/forwarded.txt"
''')

    def call(self, *options):
        return subprocess.run([str(self.repo / 'bin/MAGverify'), *map(str, options)],
                              cwd=self.repo, env=self.env, text=True, capture_output=True)

    def test_help_for_every_module(self):
        for name in STAGES:
            with self.subTest(name=name):
                result = self.call(name, '--help')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('--output-dir', result.stdout)
                self.assertIn('--contigs-dir', result.stdout)
                self.assertNotIn('--source-assembly', result.stdout)
                self.assertNotIn('--anvio-env', result.stdout)
                self.assertIn('--output-directory', result.stdout)  # legacy alias

    def test_legacy_output_directory_option(self):
        p = self.call('map', '--contigs-dir', self.source_dir, '--reads', self.reads,
                      '--output-directory', self.output, '--threads', '2')
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue((self.output / 'results' / 'mapping_manifest.tsv').is_symlink())

    def test_clipping_direct_files_without_study_directory(self):
        p = self.call('clipping', '-i', self.assembly, '--mag-fasta', self.mag,
                      '--reads', self.reads, '--contigs-dir', self.source_dir,
                      '-o', self.output, '-t', '2')
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue((self.output / '01_mapping/assembler_mapping_outputs.tsv').is_file())
        self.assertEqual((self.output / '00_prep/mag_fasta_links/Bin13L.fna').resolve(), self.mag)
        self.assertEqual((self.output / '01_mapping/assembler_mapping_outputs.tsv').is_file(), True)
        self.assertTrue((self.output / '03_clipping/forwarded.txt').is_file())

    def test_downstream_reuses_mapping_with_reads_automatically(self):
        assemblies = self.home / 'assembler_fastas'
        assemblies.mkdir()
        shutil.copy2(self.source, assemblies / self.source.name)
        first = self.call('map', '--reads', self.reads, '--contigs-dir', assemblies,
                          '--output-dir', self.output, '--threads', '2')
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        mapping_script = self.repo / 'modules/01_mapping.sh'
        mapping_script.write_text('#!/bin/bash\nexit 77\n')
        second = self.call('clipping', '--reads', self.reads,
                           '--contigs-dir', assemblies, '--mag-fasta', self.mag,
                           '--output-dir', self.output, '--threads', '2')
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertTrue((self.output / '03_clipping/forwarded.txt').exists())

    def test_map_multiple_source_fastas_in_directory(self):
        other = self.source_dir / 'assembler_two.fa'
        other.write_text('>ctg2\nAAAA\n')
        p = self.call('map', '--reads', self.reads, '--contigs-dir', self.source_dir,
                      '-o', self.output)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        rows = (self.output / '01_mapping/assembler_mapping_outputs.tsv').read_text().splitlines()
        self.assertEqual(len(rows), 3)
        self.assertIn('assembler_two.fa', rows[2])


    def test_reuses_an_explicit_manifest(self):
        p = self.call('clipping', '-i', self.assembly, '--mag-fasta', self.mag,
                      '--mapping-manifest', self.manifest, '-o', self.output)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertFalse((self.output / '01_mapping').exists())

    def test_clipping_directories_build_combined_automatically(self):
        mags = self.home / 'my_MAG_fastas'
        mags.mkdir()
        shutil.copy2(self.mag, mags / self.mag.name)
        assemblies = self.home / 'original_assemblies'
        assemblies.mkdir()
        shutil.copy2(self.source, assemblies / self.source.name)
        p = self.call('clipping', '--contigs-dir', assemblies, '--mag-dir', mags,
                      '--reads', self.reads, '--output-dir', self.output,
                      '--threads', '2')
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        combined = self.output / '00_prep/combined_assembly.fasta'
        self.assertEqual(combined.read_text(), self.source.read_text())
        self.assertTrue((self.output / '01_mapping/assembler_mapping_outputs.tsv').is_file())
        forwarded = (self.output / '03_clipping/forwarded.txt').read_text()
        self.assertIn(str(mags), forwarded)

    def test_automatic_combination_from_existing_manifest(self):
        mags = self.home / 'MAGs'
        mags.mkdir()
        shutil.copy2(self.mag, mags / self.mag.name)
        p = self.call('clipping', '--mag-dir', mags,
                      '--mapping-manifest', self.manifest,
                      '--output-dir', self.output)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual((self.output / '00_prep/combined_assembly.fasta').read_text(),
                         self.source.read_text())
        self.assertFalse((self.output / '01_mapping').exists())

    def test_reject_duplicate_contig_ids_in_assembly_directory(self):
        originals = self.home / 'originals'
        originals.mkdir()
        shutil.copy2(self.source, originals / 'first.fa')
        shutil.copy2(self.source, originals / 'second.fa')
        p = self.call('clipping', '--contigs-dir', originals, '--mag-fasta', self.mag,
                      '--reads', self.reads, '--output-dir', self.output)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('Duplicate contig ID', p.stderr)

    def test_combined_assembly_recreated_after_source_change(self):
        originals = self.home / 'orig'
        originals.mkdir()
        shutil.copy2(self.source, originals / 'source.fa')
        base = ('clipping', '--contigs-dir', originals, '--mag-fasta', self.mag,
                '--reads', self.reads, '--output-dir', self.output)
        first = self.call(*base)
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
        second = self.call(*base)
        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
        self.assertIn('(reused)', second.stdout)
        (originals / 'source.fa').write_text('>ctg1\nACGTACGT\n>ctg2\nAAAA\n')
        third = self.call(*base)
        self.assertEqual(third.returncode, 0, third.stdout + third.stderr)
        self.assertIn('>ctg2', (self.output / '00_prep/combined_assembly.fasta').read_text())

    def test_missing_source_input_explained(self):
        p = self.call('clipping', '-i', self.assembly, '--mag-fasta', self.mag,
                      '--reads', self.reads, '-o', self.output)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('--contigs-dir', p.stderr)

    def test_rejects_study_directory_as_input(self):
        p = self.call('clipping', '-i', self.data, '--mag-fasta', self.mag,
                      '-o', self.output)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('--input', p.stderr)

    def test_dry_run_does_not_write_directories(self):
        p = self.call('clipping', '-i', self.assembly, '--mag-fasta', self.mag,
                      '--reads', self.reads, '--contigs-dir', self.source_dir,
                      '-o', self.output, '--dry-run')
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn('01_mapping.sh', p.stdout)
        self.assertIn('03_clipping.sh', p.stdout)
        self.assertFalse(self.output.exists())

    def test_missing_directory_rejected(self):
        p = self.call('map', '-r', self.reads, '--contigs-dir', self.home / 'missing',
                      '-o', self.output)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('--contigs-dir', p.stderr)

    def test_old_assemblies_dir_option_is_rejected(self):
        p = self.call('map', '--reads', self.reads, '--assemblies-dir', self.source_dir,
                      '-o', self.output)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('unrecognized arguments', p.stderr)

    def test_removed_source_assembly_option_is_rejected(self):
        p = self.call('map', '--reads', self.reads, '--source-assembly', self.source,
                      '-o', self.output)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('unrecognized arguments', p.stderr)

    def test_source_dir_is_supported(self):
        p = self.call('map', '--reads', self.reads, '--contigs-dir', self.source_dir,
                      '-o', self.output, '--dry-run')
        self.assertEqual(p.returncode, 0, p.stderr + p.stdout)


    def test_each_stage_dry_run(self):
        base = ['-i', str(self.assembly), '--mag-map', str(self.map_table),
                '-o', str(self.output), '--dry-run']
        mapped = ['--mapping-manifest', str(self.manifest)]
        cases = {
            'map': ['-r', str(self.reads), '--contigs-dir', str(self.source_dir), '-o', str(self.output)],
            'clipping': ['-i', str(self.assembly), '--mag-fasta', str(self.mag),
                         '-o', str(self.output)] + mapped,
            'chimera': base + mapped + ['--taxonomy', str(self.generic_tsv),
                                        '--clipping-table', str(self.generic_tsv)],
            'circularity': ['--mag-fasta', str(self.mag), '--reads', str(self.reads),
                            '--contigs-dir', str(self.source_dir), '-o', str(self.output)],
            'false-closure': base + mapped + ['--clipping-table', str(self.generic_tsv),
                        '--clipping-all-table', str(self.generic_tsv), '--scg-table', str(self.generic_tsv),
                        '--assembler-circularity', str(self.generic_tsv)],
            'variants': base + mapped + ['--reads', str(self.reads)],
            'phantom': base + mapped + ['--reads', str(self.reads)],
            'repeats': base + ['--contigs-dir', str(self.source_dir)],
            'false-duplication': base + mapped + ['--clipping-table', str(self.generic_tsv)],
            'summary': base + ['--include', 'repeats', '--repeats-table', str(self.generic_tsv),
                                 '--repeats-mag-table', str(self.generic_tsv)],
        }
        for stage, cli in cases.items():
            with self.subTest(stage=stage):
                p = self.call(stage, *cli, '--dry-run')
                self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
                self.assertIn('.sh', p.stdout)

    def test_directory_interface_across_all_modules(self):
        source_dir = self.home / 'source_fastas'
        mag_dir = self.home / 'mag_fastas'
        source_dir.mkdir()
        mag_dir.mkdir()
        shutil.copy2(self.source, source_dir / self.source.name)
        shutil.copy2(self.mag, mag_dir / self.mag.name)
        shared = ['--contigs-dir', str(source_dir), '-o', str(self.output), '--dry-run']
        mags = ['--mag-dir', str(mag_dir)]
        manifest = ['--mapping-manifest', str(self.manifest)]
        reads = ['--reads', str(self.reads)]
        case = {
            'map': shared + reads,
            'clipping': shared + mags + reads,
            'chimera': shared + mags + manifest + ['--taxonomy', str(self.generic_tsv), '--clipping-table', str(self.generic_tsv)],
            'circularity': shared + mags + reads,
            'false-closure': shared + mags + manifest + ['--clipping-table', str(self.generic_tsv), '--clipping-all-table', str(self.generic_tsv), '--scg-table', str(self.generic_tsv), '--assembler-circularity', str(self.generic_tsv)],
            'variants': shared + mags + manifest + reads,
            'phantom': shared + mags + manifest + reads,
            'repeats': shared + mags,
            'false-duplication': shared + mags + manifest + ['--clipping-table', str(self.generic_tsv)],
            'summary': shared + mags + ['--include', 'repeats', '--repeats-table', str(self.generic_tsv), '--repeats-mag-table', str(self.generic_tsv)],
        }
        self.assertEqual(set(case), set(STAGES))
        for stage, options in case.items():
            with self.subTest(stage=stage):
                p = self.call(stage, *options)
                self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
                self.assertIn('.sh', p.stdout)

    def test_chimera_accepts_database_without_taxonomy_table(self):
        args = ('chimera', '--input', self.assembly, '--mag-map', self.map_table,
                '--mapping-manifest', self.manifest,
                '--database', self.prefix, '--clipping-table', self.generic_tsv,
                '-o', self.output, '--threads', '2', '--dry-run')
        p = self.call(*args)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn('07_chimera.sh', p.stdout)

    def test_chimera_explains_missing_taxonomy_database(self):
        p = self.call('chimera', '-i', self.assembly, '--mag-map', self.map_table,
                      '--mapping-manifest', self.manifest,
                      '--clipping-table', self.generic_tsv,
                      '-o', self.output, '--dry-run')
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('--database', p.stderr)

    def test_removed_workflow_and_preprocessing_commands(self):
        for name in ('prepare', 'structure', 'sequence', 'duplication', 'report',
                     'run', 'taxonomy', 'scg'):
            with self.subTest(name=name):
                p = self.call(name, '--help')
                self.assertNotEqual(p.returncode, 0)
                self.assertIn('invalid choice', p.stderr)

    def test_missing_stage_specific_table_explained(self):
        p = self.call('false-closure', '-i', self.assembly,
                      '--mag-map', self.map_table, '--mapping-manifest', self.manifest,
                      '-o', self.output, '--dry-run')
        self.assertNotEqual(p.returncode, 0)
        self.assertIn('--clipping-table', p.stderr)


if __name__ == '__main__':
    unittest.main()
