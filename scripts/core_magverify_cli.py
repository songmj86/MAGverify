#!/usr/bin/env python3
"""Command-line entry point for MAGverify."""
import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


class OptionHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Show options without value placeholders, each option on one line."""

    def __init__(self, prog):
        width = max(100, shutil.get_terminal_size((120, 24)).columns)
        super().__init__(prog, max_help_position=40, width=width)

    def _format_args(self, action, default_metavar):
        return '' if action.option_strings else super()._format_args(action, default_metavar)

    def _format_action_invocation(self, action):
        return ', '.join(action.option_strings) or super()._format_action_invocation(action)

    def _format_usage(self, usage, actions, groups, prefix):
        text = super()._format_usage(usage, actions, groups, prefix)
        text = re.sub(r'(?<=\S)  +(?=\S)', ' ', text.replace(' ]', ']'))
        return '\n'.join(line.rstrip() for line in text.split('\n'))


help_formatter = OptionHelpFormatter

ROOT = Path(__file__).resolve().parent.parent / "modules"
SCRIPTS_DIR = Path(__file__).resolve().parent


class _PublicTables:

    def _module(self):
        path = SCRIPTS_DIR / 'core_public_tables.py'
        if not path.is_file():
            return None
        import importlib.util
        spec = importlib.util.spec_from_file_location('magverify_public_tables', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def migrate(self, root):
        module = self._module()
        if module and not DRY_RUN:
            module.migrate(root)

    def publish(self, root):
        module = self._module()
        if module and not DRY_RUN:
            module.publish(root)


public_tables = _PublicTables()
STAGES = {
    'map': ('01_mapping.sh', '01_mapping', 'Map HiFi reads to each assembly'),
    'circularity': ('02_circularity.sh', '02_circularity', 'Validate circular genome topology'),
    'clipping': ('03_clipping.sh', '03_clipping', 'Assess read-clipping sites'),
    'false-closure': ('04_false_closure.sh', '04_false_closure', 'Assess false genome closure'),
    'phantom': ('05_phantom.sh', '05_phantom', 'Detect phantom sequences'),
    'repeats': ('06_repeats.sh', '06_repeats', 'Characterize repeat regions'),
    'chimera': ('07_chimera.sh', '07_chimera', 'Assess gene-level taxonomic chimeras'),
    'false-duplication': ('08_false_duplication.sh', '08_false_duplication', 'Assess false duplication'),
    'variants': ('09_assembly_variants.sh', '09_assembly_variants', 'Assess assembly variants'),
    'summary': ('10_summary.sh', '.cache/summary', 'Summarize and validate QC results'),
}

FILE_OPTIONS = {
    'mag_map': ('MAG_CONTIG_MAP', '--mag-map'),
    'mapping_manifest': ('ASSEMBLER_MAPPING_MANIFEST', '--mapping-manifest'),
    'taxonomy': ('GTDB_GENE_TAX_TSV', '--taxonomy'),
    'mag_taxonomy': ('MAG_GTDB_TAXONOMY_TSV', '--mag-taxonomy'),
    'clipping_table': ('CLIP_FINAL', '--clipping-table'),
    'clipping_all_table': ('CLIP_ALL_FINAL', '--clipping-all-table'),
    'scg_table': ('SCG_TSV', '--scg-table'),
    'assembler_circularity': ('CIRCULARITY_TSV', '--assembler-circularity'),
    'circularity_table': ('CIRCULARITY_ONLY_FINAL', '--circularity-table'),
    'circularity_audit': ('CIRCULARITY_ONLY_AUDIT', '--circularity-audit'),
    'false_closure_table': ('FALSE_CLOSURE_FINAL', '--false-closure-table'),
    'false_closure_audit': ('FALSE_CLOSURE_AUDIT', '--false-closure-audit'),
    'chimera_table': ('GTDB_CHIMERA_REPORT_CONTIG', '--chimera-table'),
    'variants_table': ('ALLELE_CONTIG', '--variants-table'),
    'phantom_table': ('PHANTOM_FINAL', '--phantom-table'),
    'repeats_table': ('REPEAT_SUM', '--repeats-table'),
    'repeats_mag_table': ('REPEAT_MAG_SUM', '--repeats-mag-table'),
    'duplication_table': ('FALSE_DUP_CONTIG_SUM', '--duplication-table'),
    'duplication_clipping_table': ('FALSE_DUP_CLIPPING_CONTIG_SUM', '--duplication-clipping-table'),
    'duplication_clipping_mag_table': ('FALSE_DUP_CLIPPING_MAG_SUM', '--duplication-clipping-mag-table'),
}
RESULTS = {
    'taxonomy': '00_prep/taxonomy.gene.tsv',
    'mag_map': '00_prep/MAG_contig_map.tsv',
    'mag_taxonomy': '00_prep/gtdbtk_taxonomy.MAG.tsv',
    'database': '00_prep/GTDB_R232',
    'mapping_manifest': '01_mapping/.cache/mapping_manifest.tsv',
    'clipping_table': '03_clipping/.cache/clipping.site.internal.tsv',
    'clipping_all_table': '03_clipping/.cache/clipping.site.tsv',
    'scg_table': '00_prep/scg.contig.tsv',
    'assembler_circularity': '00_prep/assembler_circularity.MAG.tsv',
    'circularity_table': '02_circularity/.cache/circularity.MAG.tsv',
    'circularity_audit': '02_circularity/.cache/circularity.MAG.audit.tsv',
    'false_closure_table': '04_false_closure/.cache/false_closure.MAG.tsv',
    'false_closure_audit': '04_false_closure/.cache/false_closure.MAG.audit.tsv',
    'chimera_table': '07_chimera/.cache/chimera.contig.tsv',
    'variants_table': '09_assembly_variants/.cache/variants.contig.tsv',
    'phantom_table': '05_phantom/.cache/phantom.region.tsv',
    'repeats_table': '06_repeats/.cache/repeats.contig.tsv',
    'repeats_mag_table': '06_repeats/.cache/repeats.MAG.tsv',
    'duplication_table': '08_false_duplication/.cache/false_duplication.contig.tsv',
    'duplication_clipping_table': '08_false_duplication/.cache/clipping/false_duplication.contig.clipping.tsv',
    'duplication_clipping_mag_table': '08_false_duplication/.cache/clipping/false_duplication.MAG.clipping.tsv',
}
FINAL_TABLES = {
    'map': ('01_mapping/mapping.assembler.tsv',),
    'circularity': ('02_circularity/circularity.MAG.tsv', '02_circularity/circularity.contig.tsv'),
    'clipping': ('03_clipping/clipping.site.tsv', '03_clipping/clipping.contig.tsv', '03_clipping/clipping.MAG.tsv',
                 '03_clipping/clipping.assembler.tsv'),
    'false-closure': ('04_false_closure/false_closure.MAG.tsv', '04_false_closure/false_closure.contig.tsv'),
    'phantom': ('05_phantom/phantom.region.tsv', '05_phantom/phantom.contig.tsv', '05_phantom/phantom.MAG.tsv'),
    'repeats': ('06_repeats/repeats.contig.tsv', '06_repeats/repeats.MAG.tsv',
                '06_repeats/repeats.region.tsv', '06_repeats/repeats.HSP.tsv'),
    'chimera': ('07_chimera/chimera.region.tsv', '07_chimera/chimera.contig.tsv', '07_chimera/chimera.MAG.tsv'),
    'false-duplication': ('08_false_duplication/false_duplication.pair.tsv',
                          '08_false_duplication/false_duplication.contig.tsv',
                          '08_false_duplication/false_duplication.MAG.tsv'),
    'variants': ('09_assembly_variants/variants.site.tsv', '09_assembly_variants/variants.contig.tsv',
                 '09_assembly_variants/variants.MAG.tsv'),
    'summary': ('summary.MAG.tsv', 'summary.contig.tsv'),
}
GENERATED_RESULTS_NAMES = {
    'mapping_manifest.tsv', 'clipping_sites.tsv', 'clipping_mag.tsv', 'clipping_contig.tsv', 'gene_taxonomy.tsv',
    'scg.tsv', 'chimera_contig.tsv', 'chimera_mag.tsv', 'circularity_mag.tsv', 'false_closure_mag.tsv',
    'false_closure_contig.tsv', 'variants_contig.tsv', 'variants_mag.tsv', 'phantom_contig.tsv', 'phantom_mag.tsv',
    'repeats_contig.tsv', 'repeats_mag.tsv', 'false_duplication_contig.tsv', 'false_duplication_mag.tsv',
    'summary_mag.tsv', 'summary_contig.tsv', 'clipping.site.tsv', 'clipping.contig.tsv', 'clipping.MAG.tsv',
    'gene_taxonomy.gene.tsv', 'scg.gene.tsv', 'chimera.contig.tsv', 'chimera.MAG.tsv', 'circularity.MAG.tsv',
    'false_closure.contig.tsv', 'false_closure.MAG.tsv', 'variants.contig.tsv', 'variants.MAG.tsv',
    'phantom.region.tsv', 'phantom.MAG.tsv', 'repeats.contig.tsv', 'repeats.MAG.tsv',
    'false_duplication.contig.tsv', 'false_duplication.MAG.tsv', 'summary.contig.tsv', 'summary.MAG.tsv',
}


def remove_generated_results_dir(root):
    """Delete OUTPUT/results/ from earlier versions when it holds only generated copies or links."""
    folder = root / 'results'
    if not folder.is_dir():
        return
    entries = list(folder.iterdir())
    if all(entry.name in GENERATED_RESULTS_NAMES and (entry.is_symlink() or entry.is_file()) for entry in entries):
        for entry in entries:
            entry.unlink()
        folder.rmdir()


def report_final_tables(stage, root):
    tables = [root / relative for relative in FINAL_TABLES.get(stage, ())]
    tables = [path for path in tables if path.is_file() and path.stat().st_size > 0]
    if tables:
        print('Final tables:', flush=True)
        for path in tables:
            print('  ' + str(path), flush=True)


FLAGS = {
    'clipping': 'RUN_CLIPPING', 'chimera': 'RUN_TAXONOMIC_CHIMERA',
    'circularity': 'RUN_CIRCULARITY', 'false-closure': 'RUN_FALSE_CLOSURE',
    'variants': 'RUN_ASSEMBLY_VARIANTS', 'phantom': 'RUN_PHANTOM',
    'repeats': 'RUN_REPEAT', 'false-duplication': 'RUN_FALSE_DUPLICATION',
}
NEEDS_MAP = {'clipping', 'chimera', 'false-closure', 'variants', 'phantom', 'false-duplication'}
NEEDS_MAG_MAP = {'chimera', 'false-closure', 'variants', 'phantom', 'repeats', 'false-duplication', 'summary'}
NEEDS_SOURCES = {'map', 'circularity', 'repeats'}

DRY_RUN = os.environ.get('MAGVERIFY_DRY_RUN') == '1'


class UsageError(Exception):
    pass


def file_input(value, name):
    if not value:
        raise UsageError(f'Missing {name}. Run MAGverify MODULE -h for input options.')
    path = Path(value).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise UsageError(f'{name} must be an existing non-empty file: {path}')
    return path


def dir_input(value, name):
    path = Path(value).expanduser().resolve()
    if path.is_file():
        raise UsageError(f'{name} needs a directory, not a file: {path} (give {path.parent})')
    if not path.is_dir():
        raise UsageError(f'{name} directory not found: {path}')
    return path


def assembly_fasta(value):
    path = file_input(value, '--input')
    if path.suffix.lower() not in ('.fa', '.fasta', '.fna'):
        raise UsageError('--input must be an uncompressed assembly FASTA (.fa, .fna, .fasta), not a study directory')
    return path


READ_SUFFIXES = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')


def reads_input(directory, root, dry_run):
    """HiFi reads from a directory. One FASTQ file is used as is; several files (all gzipped or all plain) are joined once into OUTPUT/00_prep/hifi_reads.fastq[.gz]."""
    parent = dir_input(directory, '-r')
    files = sorted(p for p in parent.iterdir()
                   if p.is_file() and p.name.lower().endswith(READ_SUFFIXES) and p.stat().st_size)
    if not files:
        raise UsageError(f'No HiFi FASTQ files ({", ".join(READ_SUFFIXES)}) found in -r {parent}')
    print(f'HiFi read files ({len(files)}): ' + ', '.join(p.name for p in files), flush=True)
    if len(files) == 1:
        return files[0]
    gzipped = {p.name.lower().endswith('.gz') for p in files}
    if len(gzipped) > 1:
        raise UsageError(f'-r {parent} mixes gzipped and uncompressed FASTQ files; use one format')
    output = root / '00_prep' / ('hifi_reads.fastq.gz' if gzipped.pop() else 'hifi_reads.fastq')
    fingerprint = [dict(path=str(p), size=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns) for p in files]
    metadata = inputs_record(output)
    if output.is_file() and metadata.is_file():
        try:
            if json.loads(metadata.read_text(encoding='utf-8')) == fingerprint:
                print(f'Joined HiFi reads: {output} (reused)', flush=True)
                return output
        except (ValueError, OSError):
            pass
    print(f'Joined HiFi reads: {output}', flush=True)
    if dry_run:
        return output
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = output.with_name('.' + output.name + '.tmp')
    with temp.open('wb') as dest:
        for source in files:
            with source.open('rb') as stream:
                shutil.copyfileobj(stream, dest, 16 * 1024 * 1024)
    os.replace(temp, output)
    metadata.parent.mkdir(parents=True, exist_ok=True)
    metadata.write_text(json.dumps(fingerprint), encoding='utf-8')
    return output


FASTA_SUFFIXES = ('.fa', '.fna', '.fasta')


def fasta_paths(directory, label):
    parent = dir_input(directory, label)
    paths = sorted(p for p in parent.iterdir() if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES)
    if not paths:
        raise UsageError(f'No FASTA files (*.fa) found in {parent}; check {label}')
    return paths


def validate_fasta_paths(paths, option):
    resolved = [file_input(p, option) for p in paths]
    if any(p.suffix.lower() not in FASTA_SUFFIXES for p in resolved):
        raise UsageError(f'{option} accepts *.fa, *.fasta and *.fna files')
    basenames = [p.name for p in resolved]
    stems = [p.stem for p in resolved]
    if len(set(basenames)) != len(paths) or len(set(stems)) != len(paths):
        raise UsageError(f'{option} FASTAs need unique basenames and IDs')
    if any(any(c.isspace() for c in p.name) for p in resolved):
        raise UsageError(f'{option} FASTA filenames cannot include whitespace')
    return resolved


def link_inputs(named_paths, directory, dry_run):
    """Make stable links within the output root; users never organize or copy inputs."""
    for name, source in named_paths:
        link = directory / name
        if link.is_symlink() and link.resolve() == source:
            continue
        if link.exists() or link.is_symlink():
            raise UsageError(f'Conflicting managed input link: {link}. Use a different output directory or remove this link.')
        if not dry_run:
            directory.mkdir(parents=True, exist_ok=True)
            link.symlink_to(source)
    return directory


def source_inputs(args, output_root):
    """Accept the directory of original, unbinned assembler FASTAs."""
    if not getattr(args, 'source_assemblies', None):
        return None
    paths = validate_fasta_paths(
        fasta_paths(args.source_assemblies, '-a'), '-a'
    )
    parent = dir_input(args.source_assemblies, '-a')
    names = [p.name for p in paths]
    print(f'Assemblies ({len(names)}): ' + ', '.join(names), flush=True)
    return parent, ' '.join(names)


def mag_inputs(args, output_root):
    if getattr(args, 'mag_dir', None):
        paths = validate_fasta_paths(fasta_paths(args.mag_dir, '-g'), '-g')
        parent = dir_input(args.mag_dir, '-g')
    else:
        return None
    print(f'MAG FASTAs ({len(paths)}): ' + ', '.join(path.name for path in paths), flush=True)
    return parent


def run(cmd, env, dry_run):
    print('$ ' + shlex.join(map(str, cmd)), flush=True)
    if not dry_run:
        start = time.monotonic()
        subprocess.run(cmd, cwd=ROOT, env=env, check=True)
        print('Elapsed: ' + format_duration(time.monotonic() - start), flush=True)


def format_duration(seconds):
    seconds = int(seconds)
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    return f'{hours}h{minutes:02d}m' if hours else f'{minutes}m{secs:02d}s'


def check_programs(stage, mapping=False, env=None):
    core = {'bash', 'python3'}
    if stage in ('map', 'clipping'):
        core |= {'samtools', 'sha256sum'}
    if stage == 'map' or mapping:
        core.add('minimap2')
    if stage == 'clipping':
        core.add('anvi-script-find-misassemblies')
    if stage == 'prep-chimera':
        core.update({'mmseqs', 'awk', 'sed'})
    if stage in ('prep-scg', 'chimera', 'variants'):
        core.add('pyrodigal')
    if stage == 'prep-scg':
        core.update({'anvi-script-reformat-fasta', 'anvi-gen-contigs-database', 'anvi-run-hmms'})
    if stage == 'variants':
        core.update({'anvi-profile', 'anvi-script-reformat-fasta',
                     'anvi-script-reformat-bam', 'anvi-gen-contigs-database',
                     'anvi-export-table'})
    extra = []
    if env:
        for var in ('ANVIO_ENV_PREFIX', 'MINIMAP2_ENV', 'MERYL_ENV'):
            if env.get(var):
                extra.append(Path(env[var]) / 'bin')
    missing = sorted(p for p in core if not shutil.which(p) and
                     not any((directory / p).is_file() for directory in extra))
    if missing:
        raise UsageError('Missing programs: ' + ', '.join(missing) + '. See the Installation section of README.md')


def check_manifest(value):
    path = file_input(value, '--mapping-manifest')
    with path.open(encoding='utf-8-sig', newline='') as handle:
        rows = csv.DictReader(handle, delimiter='\t')
        necessary = {'Assembler', 'Assembly_FASTA', 'HiFi_all_BAM'}
        if not necessary.issubset(rows.fieldnames or []):
            raise UsageError(f'--mapping-manifest lacks required columns: {sorted(necessary)}')
        count = 0
        for row in rows:
            file_input(row['Assembly_FASTA'], f'Assembly_FASTA in {path}')
            bam = file_input(row['HiFi_all_BAM'], f'HiFi_all_BAM in {path}')
            file_input(str(Path(bam).parent / '.cache' / (Path(bam).name + '.inputs.sha256')), f'BAM mapping provenance in {path}')
            count += 1
        if not count:
            raise UsageError(f'--mapping-manifest is empty: {path}')
    return path


def inputs_record(path):
    """Sidecar JSON recording the inputs of a generated file."""
    path = Path(path)
    record = path.parent / '.cache' / (path.name + '.inputs.json')
    for legacy in (Path(str(path) + '.inputs.json'), Path(str(path) + '.sources.json')):
        if legacy.is_file() and not record.exists():
            if DRY_RUN:
                return legacy
            record.parent.mkdir(exist_ok=True)
            os.replace(legacy, record)
    return record


def all_assemblies_path(root):
    """Path of the concatenated assembler FASTAs."""
    prep = root / '00_prep'
    output = prep / 'all_assemblies.fa'
    for legacy in (prep / 'all_assemblies.fasta', prep / 'combined_assembly.fasta'):
        if not legacy.is_file() or output.exists():
            continue
        if DRY_RUN:
            return legacy
        for suffix in ('', '.fai', '.sources.json', '.inputs.json'):
            old = Path(str(legacy) + suffix)
            if old.exists():
                os.replace(old, Path(str(output) + suffix))
        old_record = prep / '.cache' / (legacy.name + '.inputs.json')
        if old_record.is_file():
            os.replace(old_record, prep / '.cache' / (output.name + '.inputs.json'))
        inputs_record(output)
    return output


def assemble_combined(paths, output, dry_run):
    """Build a FASTA spanning the original assemblies, retaining contig IDs."""
    sources = [file_input(p, 'assembler FASTA') for p in paths]
    if not sources:
        raise UsageError('No assembler FASTAs available to build the all-assemblies FASTA')
    fingerprint = [dict(path=str(p.resolve()), size=p.stat().st_size,
                        mtime_ns=p.stat().st_mtime_ns) for p in sources]
    metadata = inputs_record(output)
    if output.is_file() and output.stat().st_size and metadata.is_file():
        try:
            if json.loads(metadata.read_text(encoding='utf-8')) == fingerprint:
                print(f'All assemblies: {output} (reused)', flush=True)
                return output
        except (ValueError, OSError):
            pass
    print(f'All assemblies: {output} ({len(sources)} assembler FASTAs)', flush=True)
    if dry_run:
        return output
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix='.combined_', suffix='.fa', dir=output.parent)
    temp = Path(temp_name)
    identifiers = set()
    try:
        with os.fdopen(fd, 'wb') as dest:
            for source in sources:
                seen_header = False
                with source.open('rb') as stream:
                    for line in stream:
                        if line.startswith(b'>'):
                            fields = line[1:].split(maxsplit=1)
                            if not fields:
                                raise UsageError(f'Empty FASTA header in {source}')
                            try:
                                identifier = fields[0].decode('utf-8')
                            except UnicodeDecodeError:
                                raise UsageError(f'Invalid FASTA contig identifier in {source}') from None
                            if identifier in identifiers:
                                raise UsageError(f'Duplicate contig ID {identifier!r} across assemblies. '
                                                 'IDs must be unique; inspect your assembler FASTAs.')
                            identifiers.add(identifier)
                            seen_header = True
                        elif line.strip() and not seen_header:
                            raise UsageError(f'Sequence data appears before a FASTA header in {source}')
                        dest.write(line if line.endswith(b'\n') else line + b'\n')
                if not seen_header:
                    raise UsageError(f'No FASTA records in {source}')
        os.replace(temp, output)
        metadata.parent.mkdir(parents=True, exist_ok=True)
        metadata.write_text(json.dumps(fingerprint, indent=2) + '\n', encoding='utf-8')
        Path(str(output) + '.fai').unlink(missing_ok=True)
    finally:
        temp.unlink(missing_ok=True)
    return output


def source_fastas_from_manifest(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        rows = csv.DictReader(stream, delimiter='\t')
        if not rows.fieldnames or 'Assembly_FASTA' not in rows.fieldnames:
            raise UsageError(f'Mapping manifest lacks Assembly_FASTA: {path}')
        paths = [file_input(row['Assembly_FASTA'], 'Assembly_FASTA in mapping manifest') for row in rows]
    if not paths:
        raise UsageError(f'No assembler FASTAs in mapping manifest: {path}')
    if len(set(paths)) != len(paths):
        raise UsageError('Mapping manifest lists the same assembler FASTA more than once')
    return paths


def source_paths(sources):
    directory, basenames = sources
    return [directory / name for name in basenames.split()]


def previously_mapped_sources(args, root):
    """Use original assembler references recorded by a previous mapping run."""
    manifest = getattr(args, 'mapping_manifest', None)
    if manifest:
        path = check_manifest(manifest)
    else:
        path = root / RESULTS['mapping_manifest']
        if not path.is_file():
            return None
        path = check_manifest(path)
    return source_fastas_from_manifest(path)


def source_links_from_manifest(args, root):
    paths = previously_mapped_sources(args, root)
    if not paths:
        return None
    named = []
    used = set()
    for index, path in enumerate(paths, start=1):
        name = path.name
        if name in used:
            name = f'assembler_{index}_{name}'
        used.add(name)
        named.append((name, path))
    for legacy in (root / '00_prep' / 'source_assembly_links', root / '00_prep' / 'assembler_links'):
        if (not args.dry_run and legacy.is_dir()
                and all(p.is_symlink() or p.name.endswith('.fai') for p in legacy.iterdir())):
            shutil.rmtree(legacy)
    directory = link_inputs(named, root / '00_prep' / '.cache' / 'assembler_links', args.dry_run)
    return directory, ' '.join(name for name, _ in named)


OUTPUT_FILE_SUFFIXES = ('.tsv', '.txt', '.csv', '.fa', '.fasta', '.fna', '.gz', '.bam')


def defaults_for(stage, args):
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    threads = getattr(args, 'threads', None)
    threads = threads if threads is not None else min(100, cpus)
    if threads < 1:
        raise UsageError('-t must be a positive integer')
    if threads > cpus:
        raise UsageError(f'-t {threads} exceeds the available CPUs ({cpus})')
    root = Path(args.output_directory).expanduser().resolve()
    if root.name in {'00_prep', '10_summary'} | {folder for _script, folder, _text in STAGES.values()}:
        print(f'-o is the output folder; using {root.parent} (given: {root})', flush=True)
        root = root.parent
    if root.is_file() or root.suffix.lower() in OUTPUT_FILE_SUFFIXES:
        raise UsageError(f'-o is the output folder shared by all commands, not a file: {root}. '
                         'Each command writes its tables inside it '
                         '(prep-taxonomy: -o DIR writes DIR/00_prep/gtdbtk_taxonomy.MAG.tsv)')
    env = os.environ.copy()
    for key in list(env):
        if key in {'OUT', 'READS', 'ASSEMBLY_FA', 'MAG_DIR', 'MAG_CONTIG_MAP',
                   'ASSEMBLER_ASSEMBLY_DIR', 'ASSEMBLER_FASTAS', 'ASSEMBLER_MAPPING_MANIFEST',
                   'MAPPING_DIR', 'CLIPPING_DIR', 'PREP_DIR', 'PER_ASSEMBLER_MANIFEST',
                   'PER_ASSEMBLER_MAP_ROOT', 'HIFI_MAPPING_DIR', 'SCG_TSV', 'SCG_OUT_DIR', 'GENE_CALLS_DIR',
                   'CIRCULARITY_TSV', 'MAG_GTDB_TAXONOMY_TSV', 'GTDB_GENE_TAX_TSV',
                   'CLIP_FINAL', 'CLIP_ALL_FINAL', 'GTDB_CHIMERA_REPORT_CONTIG',
                   'CIRCULARITY_ONLY_FINAL', 'CIRCULARITY_ONLY_AUDIT', 'FALSE_CLOSURE_FINAL',
                   'FALSE_CLOSURE_AUDIT', 'ALLELE_CONTIG', 'PHANTOM_FINAL', 'REPEAT_SUM',
                   'REPEAT_MAG_SUM', 'FALSE_DUP_CONTIG_SUM', 'FALSE_DUP_CLIPPING_CONTIG_SUM',
                   'FALSE_DUP_CLIPPING_MAG_SUM', 'MEMORY_GB'} or key.startswith((
                       'MAPPING_INPUT_', 'CLIPPING_INPUT_', 'CHIMERA_INPUT_', 'CIRCULARITY_INPUT_',
                       'FALSE_CLOSURE_INPUT_', 'VARIANT_INPUT_', 'PHANTOM_INPUT_', 'REPEAT_INPUT_',
                       'FALSE_DUP_INPUT_', 'SUMMARY_INPUT_')):
            env.pop(key, None)
    env.update(OUT=str(root), THREADS=str(threads), QC_SKIP_LOCAL_CONFIG='1')
    migrate_prep_layout(root)
    return root, env


def prepare_map(args, root, env, main_input):
    map_file = root / RESULTS['mag_map']
    if getattr(args, 'mag_map', None):
        return file_input(args.mag_map, '--mag-map')
    if 'MAG_DIR' not in env:
        if map_file.is_file():
            return file_input(map_file, 'MAG map in output directory')
        raise UsageError('MAG membership is needed: pass -g DIR or --mag-map TSV')
    if not main_input:
        raise UsageError('Preparing a MAG map needs -a DIR')
    if args.dry_run:
        run(['samtools', 'faidx', str(main_input)], env, True)
        run(['python3', str(SCRIPTS_DIR / 'prep_build_mag_contig_map.py'), '--mag-dir', env['MAG_DIR'],
             '--extensions', 'fa fasta fna', '--assembly-fai', str(main_input) + '.fai',
             '--output', str(map_file)], env, True)
        return map_file
    if not shutil.which('samtools'):
        raise UsageError('samtools is necessary to index the assembly and prepare MAG membership')
    fai = Path(str(main_input) + '.fai')
    if not fai.is_file() or fai.stat().st_mtime_ns < Path(main_input).stat().st_mtime_ns:
        run(['samtools', 'faidx', str(main_input)], env, False)
    map_file.parent.mkdir(parents=True, exist_ok=True)
    candidate = map_file.with_name(map_file.name + '.new')
    run(['python3', str(SCRIPTS_DIR / 'prep_build_mag_contig_map.py'), '--mag-dir', env['MAG_DIR'],
         '--extensions', 'fa fasta fna', '--assembly-fai', str(fai),
         '--output', str(candidate)], env, False)
    if map_file.is_file() and map_file.read_bytes() == candidate.read_bytes():
        candidate.unlink()
    else:
        os.replace(candidate, map_file)
    return file_input(map_file, 'generated MAG map')


def mapping_result(stage, args, root, env, sources):
    if args.mapping_manifest:
        return check_manifest(args.mapping_manifest)
    generated = root / RESULTS['mapping_manifest']
    if generated.is_file():
        manifest = check_manifest(generated)
        if sources:
            requested = {p.resolve() for p in source_paths(sources)}
            recorded = {p.resolve() for p in source_fastas_from_manifest(manifest)}
            if requested != recorded:
                raise UsageError('Assemblies differ from the existing mapping manifest. '
                                 'Rerun MAGverify map with the new inputs or use another -o DIR.')
        return manifest
    if sources and 'READS' not in env:
        raise UsageError('Run MAGverify map with -r first, or supply -r DIR to map automatically')
    if 'READS' in env and sources:
        if not args.dry_run:
            check_programs('map', env=env)
        run(['bash', str(ROOT / '01_mapping.sh'), env['READS'], str(root / '01_mapping')], env, args.dry_run)
        return generated if args.dry_run else check_manifest(generated)
    if getattr(args, 'reads', None) and not sources and not generated.is_file():
        raise UsageError('Provide -a DIR with -r DIR, or --mapping-manifest')
    raise UsageError(f'{stage} needs the map result in this -o DIR: run MAGverify map first')


def reads_from_manifest(args, root):
    """HiFi reads recorded by MAGverify map (Reads column), if available."""
    manifest = getattr(args, 'mapping_manifest', None)
    path = Path(manifest).expanduser() if manifest else root / RESULTS['mapping_manifest']
    if not path.is_file():
        return None
    with path.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    values = {(row.get('Reads') or '').strip() for row in rows} - {''}
    if not values:
        return None
    if len(values) > 1:
        raise UsageError(f'Mapping manifest lists more than one read file: {path}')
    reads = file_input(values.pop(), f'Reads in mapping manifest {path}')
    print(f'HiFi reads: {reads} (from mapping manifest)', flush=True)
    return reads


def choose_inputs(args, root, env):
    stage = args.module
    main_input = None
    if getattr(args, 'input', None):
        main_input = assembly_fasta(args.input)
        env['ASSEMBLY_FA'] = str(main_input)
    if getattr(args, 'reads', None):
        env['READS'] = str(reads_input(args.reads, root, args.dry_run))
    elif stage != 'map':
        recorded = reads_from_manifest(args, root)
        if recorded:
            env['READS'] = str(recorded)
    sources = source_inputs(args, root)
    if not sources and stage != 'map':
        sources = source_links_from_manifest(args, root)
    if sources:
        env['ASSEMBLER_ASSEMBLY_DIR'], env['ASSEMBLER_FASTAS'] = str(sources[0]), sources[1]
    if stage not in ('map', 'circularity') and not main_input:
        if sources:
            main_input = assemble_combined(source_paths(sources),
                                           all_assemblies_path(root), args.dry_run)
        else:
            stored = all_assemblies_path(root)
            if stored.is_file() and stored.stat().st_size:
                main_input = stored
                print(f'All assemblies: {stored} (existing)', flush=True)
        if main_input:
            env['ASSEMBLY_FA'] = str(main_input)
    mags = mag_inputs(args, root)
    if mags:
        env['MAG_DIR'] = str(mags)
    for option, (variable, flag) in FILE_OPTIONS.items():
        if option in ('mapping_manifest', 'mag_map'):
            continue
        value = getattr(args, option, None)
        if value:
            env[variable] = str(file_input(value, flag))
    conda_prefix = env.get('CONDA_PREFIX')
    if conda_prefix and (Path(conda_prefix) / 'bin' / 'anvi-profile').is_file():
        env['ANVIO_ENV_PREFIX'] = conda_prefix
    elif shutil.which('anvi-profile'):
        env['ANVIO_ENV_PREFIX'] = str(Path(shutil.which('anvi-profile')).resolve().parent.parent)
    else:
        env.pop('ANVIO_ENV_PREFIX', None)
    if stage == 'variants':
        env['ALLELE_ANVIO_FULL_PROFILE_THREADS'] = env['THREADS']
        if not env.get('ANVIO_ENV_PREFIX'):
            binary = shutil.which('anvi-profile')
            if binary:
                env['ANVIO_ENV_PREFIX'] = str(Path(binary).resolve().parent.parent)
    memory_gb = getattr(args, 'memory_gb', None)
    if memory_gb is not None and memory_gb < 1:
        raise UsageError('-m must be a positive integer')
    if memory_gb:
        env['MEMORY_GB'] = str(memory_gb)
        if stage == 'variants':
            env['ASSEMBLY_VARIANT_MEMORY_BUDGET_GB'] = str(memory_gb)
    if stage == 'variants' and args.memory_reserve_gb is not None:
        if args.memory_reserve_gb < 1:
            raise UsageError('--memory-reserve-gb must be a positive integer')
        env['ASSEMBLY_VARIANT_MEMORY_RESERVE_GB'] = str(args.memory_reserve_gb)
    if getattr(args, 'include', None):
        if stage not in ('summary'):
            raise UsageError('--include is only valid for summary')
        selected = {name.strip().lower() for name in args.include.split(',') if name.strip()}
        bad = selected - set(FLAGS)
        if bad or not selected:
            raise UsageError('--include must list one or more of: ' + ', '.join(FLAGS))
        for name, var in FLAGS.items():
            env[var] = 'Y' if name in selected else 'N'
    return main_input, sources


def need_result(name, args, root, env, stage, optional=False):
    var, flag = FILE_OPTIONS[name]
    if var in env:
        return file_input(env[var], flag)
    if name not in RESULTS:
        if optional:
            return None
        raise UsageError(f'{stage} needs {flag} FILE')
    fallback = root / RESULTS[name]
    if fallback.is_file() and fallback.stat().st_size:
        env[var] = str(fallback)
        return fallback
    if optional:
        return None
    hint = {'assembler_circularity': ' (run MAGverify prep-circularity first)',
            'mag_taxonomy': ' (run MAGverify prep-taxonomy first)',
            'clipping_table': ' (run MAGverify clipping first)',
            'clipping_all_table': ' (run MAGverify clipping first)'}.get(name, '')
    raise UsageError(f'{stage} needs {flag} FILE or an existing result at {fallback}{hint}')


def prep_chimera_run(args):
    """Build an MMseqs2 taxonomy DB from user-downloaded GTDB R232 inputs."""
    root, env = defaults_for('prep-chimera', args)
    script = ROOT / '00_prep_chimera.sh'
    if not script.is_file():
        raise UsageError(f'Missing installed GTDB preparation module: {script}')
    gtdb_dir = dir_input(args.gtdb_dir, '-d')
    inputs = ('gtdb_proteins_aa_reps_r232.tar.gz',
              'bac120_taxonomy_r232.tsv', 'ar53_taxonomy_r232.tsv')
    for name in inputs:
        file_input(gtdb_dir / name, f'GTDB R232 input in -d ({name})')
    if not args.dry_run:
        check_programs('prep-chimera', env=env)
    work = root / '00_prep' / '.cache' / 'prep_chimera_work'
    db = root / RESULTS['database']
    if db.is_dir():
        raise UsageError(f'{db} is a directory; remove it so the MMseqs2 database files can be written')
    print('GTDB R232 input directory: ' + str(gtdb_dir), flush=True)
    print('GTDB R232 database prefix: ' + str(db), flush=True)
    print('Temporary MMseqs2 directory: ' + str(work), flush=True)
    run(['bash', str(script), str(work), env['THREADS'], str(gtdb_dir), str(db)],
        env, args.dry_run)


def prep_circularity_run(args):
    """Build the assembler-reported circularity table from native assembler outputs."""
    root, env = defaults_for('prep-circularity', args)
    mags = dir_input(args.mag_dir, '-g')
    files = [str(file_input(path, '-c')) for path in args.sources]
    run(['bash', str(ROOT / '00_prep_circularity.sh'), str(mags), str(root), *files], env, args.dry_run)


def prep_taxonomy_run(args):
    """Convert GTDB-Tk summaries to the MAG taxonomy table read by chimera."""
    root, env = defaults_for('prep-taxonomy', args)
    gtdbtk = dir_input(args.gtdbtk_dir, '-d')
    mags = dir_input(args.mag_dir, '-g')
    output = root / RESULTS['mag_taxonomy']
    print('MAG taxonomy: ' + str(output), flush=True)
    run(['bash', str(ROOT / '00_prep_taxonomy.sh'), str(gtdbtk), str(mags), str(output)], env, args.dry_run)


def prep_scg_run(args):
    """Annotate single-copy genes on all assemblies with anvi'o."""
    root, env = defaults_for('prep-scg', args)
    sources = source_inputs(args, root)
    if not sources:
        raise UsageError('prep-scg needs -a DIR')
    assembly = assemble_combined(source_paths(sources),
                                 all_assemblies_path(root), args.dry_run)
    if not args.dry_run:
        check_programs('prep-scg', env=env)
    out = root / '00_prep' / '.cache' / 'scg'
    env['SCG_TSV'] = str(root / RESULTS['scg_table'])
    env['GENE_CALLS_DIR'] = str(root / '00_prep')
    print('SCG table: ' + env['SCG_TSV'], flush=True)
    run(['bash', str(ROOT / '00_prep_scg.sh'), str(assembly), str(out), env['THREADS']],
        env, args.dry_run)


def add_prep_parser(subparsers):
    """A setup command, separate from the ten analysis stages."""
    p = subparsers.add_parser('prep-chimera', formatter_class=help_formatter,
        help='Build the GTDB R232 MMseqs2 database',
        description='Build the GTDB R232 MMseqs2 database', add_help=False)
    p.add_argument('-h', action='help', help=argparse.SUPPRESS)
    p.add_argument('-d', dest='gtdb_dir', required=True, metavar='DIR',
                   help='GTDB R232 directory (gtdb_proteins_aa_reps_r232.tar.gz, '
                        'bac120_taxonomy_r232.tsv, ar53_taxonomy_r232.tsv)')
    p.add_argument('-t', dest='threads', type=int, metavar='N', help='Number of threads')
    p.add_argument('-o', dest='output_directory', required=True, metavar='DIR',
                   help='Directory of base output, the same -o as the analysis modules '
                        '(database: 00_prep/GTDB_R232)')
    p = subparsers.add_parser('prep-circularity', formatter_class=help_formatter,
        help='Build the assembler circularity table',
        description='Build the assembler-reported circularity table (used by circularity and false-closure)',
        add_help=False)
    p.add_argument('-h', action='help', help=argparse.SUPPRESS)
    p.add_argument('-c', dest='sources', action='append', required=True, metavar='FILE',
                   help='Assembler circularity file; repeat for each assembler '
                        '(hifiasm-meta p_ctg.gfa, metaMDBG contigs.fasta.gz, metaFlye assembly_info.txt)')
    p.add_argument('-g', dest='mag_dir', required=True, metavar='DIR', help='Directory of MAGs (*.fa)')
    p.add_argument('-o', dest='output_directory', required=True, metavar='DIR',
                   help='Directory of base output, the same -o as the analysis modules '
                        '(table: 00_prep/assembler_circularity.MAG.tsv)')
    p = subparsers.add_parser('prep-taxonomy', formatter_class=help_formatter,
        help='Build the MAG taxonomy table from GTDB-Tk results',
        description='Build the MAG taxonomy table from GTDB-Tk results (used by chimera)',
        add_help=False)
    p.add_argument('-h', action='help', help=argparse.SUPPRESS)
    p.add_argument('-d', dest='gtdbtk_dir', required=True, metavar='DIR',
                   help='Directory of GTDB-Tk classify_wf results '
                        '(gtdbtk.bac120.summary.tsv; gtdbtk.ar53.summary.tsv)')
    p.add_argument('-g', dest='mag_dir', required=True, metavar='DIR', help='Directory of MAGs (*.fa)')
    p.add_argument('-o', dest='output_directory', required=True, metavar='DIR',
                   help='Directory of base output, the same -o as the analysis modules '
                        '(table: 00_prep/gtdbtk_taxonomy.MAG.tsv)')
    p = subparsers.add_parser('prep-scg', formatter_class=help_formatter,
        help="Annotate single-copy genes with anvi'o",
        description="Annotate single-copy genes on all assemblies with anvi'o "
                    "(used by false-closure and chimera)", add_help=False)
    p.add_argument('-h', action='help', help=argparse.SUPPRESS)
    p.add_argument('-a', dest='source_assemblies', required=True, metavar='DIR',
                   help='Directory of assembled contigs (*.fa)')
    p.add_argument('-t', dest='threads', type=int, metavar='N', help='Number of threads')
    p.add_argument('-o', dest='output_directory', required=True, metavar='DIR',
                   help='Directory of base output, the same -o as the analysis modules '
                        '(table: 00_prep/scg.contig.tsv)')
    return p

LEGACY_OUTPUT_NAMES = (
    ('MAG_phantom_sequences.final.tsv', 'phantom.MAG_region.tsv'),
    ('MAG_phantom_summary.final.tsv', 'phantom.MAG.tsv'),
    ('phantom_sequences.final.tsv', 'phantom.region.tsv'),
    ('contig_clipping.final.tsv', 'clipping.contig.tsv'),
    ('MAG_clipping.final.tsv', 'clipping.MAG.tsv'),
    ('MAG_circularity.final.tsv', 'circularity.MAG.tsv'),
    ('MAG_circularity.audit.tsv', 'circularity.MAG.audit.tsv'),
    ('MAG_END_START_closure.final.tsv', 'circularity.MAG.junction.tsv'),
    ('MAG_contig_false_closure.final.tsv', 'false_closure.MAG_contig.tsv'),
    ('contig_false_closure.final.tsv', 'false_closure.contig.tsv'),
    ('MAG_false_closure.final.tsv', 'false_closure.MAG.tsv'),
    ('MAG_false_closure.audit.tsv', 'false_closure.MAG.audit.tsv'),
    ('contig_assembly_allele_support.final.tsv', 'variants.contig.tsv'),
    ('assembly_allele_support.by_MAG.tsv', 'variants.MAG.tsv'),
    ('contig_repeat_breadth.final.tsv', 'repeats.contig.tsv'),
    ('MAG_repeat_breadth.final.tsv', 'repeats.MAG.tsv'),
    ('false_duplication.by_contig.clipping.tsv', 'false_duplication.contig.clipping.tsv'),
    ('false_duplication.by_MAG.clipping.tsv', 'false_duplication.MAG.clipping.tsv'),
    ('false_duplication.by_contig.tsv', 'false_duplication.contig.tsv'),
    ('false_duplication.by_MAG.tsv', 'false_duplication.MAG.tsv'),
    ('integrated_QC_summary.MAG_contigs.tsv', 'summary.contig.tsv'),
    ('integrated_QC_summary.by_MAG.tsv', 'summary.MAG.tsv'),
    ('contig_taxonomic_chimera.GTDB_R232.final.tsv', 'chimera.contig.tsv'),
    ('MAG_contig_taxonomic_chimera.GTDB_R232.final.tsv', 'chimera.MAG_contig.tsv'),
    ('MAG_taxonomic_chimera.GTDB_R232.final.tsv', 'chimera.MAG.tsv'),
    ('source_group_outputs.tsv', 'assembler_group_outputs.tsv'),
)
LEGACY_OUTPUT_DIRS = (
    ('02_circularity/00_source_plan', '02_circularity/00_assembler_plan'),
    ('02_circularity/01_by_source', '02_circularity/01_by_assembler'),
)

PREP_LEGACY_FILES = (
    ('genes/genes.faa', 'genes.faa'),
    ('genes/genes.gff', 'genes.gff'),
    ('genes/genes.coordinates.tsv', 'genes.coordinates.tsv'),
    ('genes/.cache/genes.inputs.sha256', '.cache/genes.inputs.sha256'),
    ('SCG/SCG_from_anvio.original_contig.tsv', 'scg.contig.tsv'),
    ('SCG/00_SCG.log', 'SCG/prep_scg.log'),
    ('gtdb_taxonomy.MAG.tsv', 'gtdbtk_taxonomy.MAG.tsv'),
    ('taxonomy.MAG.tsv', 'gtdbtk_taxonomy.MAG.tsv'),
    ('gene_taxonomy/GTDB_R232/gene_taxonomy.mmseqs2.GTDB_R232.all_contigs.tsv', 'taxonomy.gene.tsv'),
)
PREP_LEGACY_DIRS = (
    ('SCG', '.cache/scg'),
    ('gene_taxonomy/GTDB_R232', '.cache/gene_taxonomy'),
    ('tmp_prep_chimera', '.cache/prep_chimera_work'),
)


def migrate_prep_layout(root):
    """Move 00_prep files of earlier versions (genes/, SCG/, gene_taxonomy/) to the flat layout."""
    prep = root / '00_prep'
    if DRY_RUN or not prep.is_dir():
        return
    for old, new in PREP_LEGACY_FILES:
        if (prep / old).is_file() and not (prep / new).exists():
            (prep / new).parent.mkdir(parents=True, exist_ok=True)
            os.replace(prep / old, prep / new)
    for old, new in PREP_LEGACY_DIRS:
        if (prep / old).is_dir() and not (prep / new).exists():
            (prep / new).parent.mkdir(parents=True, exist_ok=True)
            os.replace(prep / old, prep / new)
    all_assemblies_path(root)
    links = prep / 'assembler_links'
    if links.is_dir() and all(p.is_symlink() or p.name.endswith('.fai') for p in links.iterdir()):
        shutil.rmtree(links)
    for leftover in ('genes', 'gene_taxonomy'):
        folder = prep / leftover
        if folder.is_dir() and not any(p.is_file() and p.name not in ('genes.lock',) for p in folder.rglob('*')):
            shutil.rmtree(folder)


def migrate_legacy_outputs(root):
    """Rename result files and folders written by earlier versions so they are reused."""
    for old, new in LEGACY_OUTPUT_DIRS:
        if (root / old).is_dir() and not (root / new).exists():
            os.replace(root / old, root / new)
    names = dict(LEGACY_OUTPUT_NAMES)
    for directory, _, files in os.walk(root):
        if Path(directory).name == 'results':
            continue
        for name in files:
            base = name[:-len('.inputs.sha256')] if name.endswith('.inputs.sha256') else name
            if base in names:
                old = Path(directory) / name
                new = Path(directory) / (names[base] + name[len(base):])
                if not new.exists():
                    os.replace(old, new)
    for module_dir in sorted(root.glob('[0-1][0-9]_*')):
        if not module_dir.is_dir() or module_dir.name == '01_mapping':
            continue
        for meta in module_dir.glob('*.inputs.sha256'):
            name = 'variants.inputs.sha256' if meta.name == 'assembly_allele_support.inputs.sha256' else meta.name
            target = module_dir / '.cache' / name
            if not target.exists():
                target.parent.mkdir(exist_ok=True)
                os.replace(meta, target)


def stage_run(args):
    stage = args.module
    root, env = defaults_for(stage, args)
    if not args.dry_run and root.is_dir():
        migrate_legacy_outputs(root)
        remove_generated_results_dir(root)
        public_tables.migrate(root)
    main_input, sources = choose_inputs(args, root, env)
    if stage not in ('map', 'circularity') and not main_input:
        raise UsageError(f'{stage} needs the map result in this -o DIR: run MAGverify map first '
                         '(or pass -a DIR)')
    if stage in NEEDS_SOURCES and not sources:
        raise UsageError(f'{stage} needs the map result in this -o DIR: run MAGverify map first '
                         '(or pass -a DIR)')
    if stage == 'map' and not args.reads:
        raise UsageError('map needs -r DIR (HiFi FASTQ/FASTQ.GZ files)')
    if stage == 'circularity' and 'MAG_DIR' not in env:
        raise UsageError('circularity needs -g DIR')
    if stage in ('circularity', 'variants', 'phantom') and 'READS' not in env:
        raise UsageError(f'{stage} needs HiFi reads: run MAGverify map first in this -o DIR, or pass -r DIR')
    if stage == 'clipping' and 'MAG_DIR' not in env:
        raise UsageError('clipping needs -g DIR (one FASTA per MAG)')
    if stage == 'chimera':
        existing_taxonomy = need_result('taxonomy', args, root, env, stage, optional=True)
        if not existing_taxonomy:
            prefix = (Path(args.database).expanduser().resolve() if args.database
                      else root / RESULTS['database'])
            if not args.database and not Path(str(prefix) + '.dbtype').is_file():
                raise UsageError(f'chimera needs the MMseqs2 database at {prefix} '
                                 '(run MAGverify prep-chimera with the same -o first), '
                                 'or --database PREFIX, or --taxonomy TSV')
            print(f'MMseqs2 database: {prefix}', flush=True)
            if not Path(str(prefix) + '.dbtype').is_file():
                raise UsageError(f'--database must be a MMseqs2 database prefix with .dbtype: {prefix}')
            env['MMSEQS_TAX_DB'] = str(prefix)
            env['GTDB_GENE_TAX_TSV'] = str(root / RESULTS['taxonomy'])
        need_result('clipping_table', args, root, env, stage)
        table = need_result('mag_taxonomy', args, root, env, stage, optional=True)
        print(f'MAG taxonomy: {table}' if table else
              'MAG taxonomy: none (SCG redundancy annotation is skipped; run MAGverify prep-taxonomy)',
              flush=True)
    if stage == 'false-closure':
        for table in ('clipping_table', 'clipping_all_table', 'assembler_circularity'):
            need_result(table, args, root, env, stage)
    if stage == 'circularity':
        table = need_result('assembler_circularity', args, root, env, stage, optional=True)
        print(f'Assembler circularity: {table}' if table else
              'Assembler circularity: none (Circular = read-supported call)', flush=True)
    if stage == 'false-duplication':
        need_result('clipping_table', args, root, env, stage)
    if getattr(args, 'mag_map', None) and stage == 'clipping':
        env['MAG_CONTIG_MAP'] = str(file_input(args.mag_map, '--mag-map'))
    if stage in NEEDS_MAG_MAP:
        env['MAG_CONTIG_MAP'] = str(prepare_map(args, root, env, main_input))
    if not args.dry_run:
        check_programs(stage, env=env)
    if stage == 'map':
        cmd = ['bash', str(ROOT / '01_mapping.sh'), env['READS'], str(root / '01_mapping')]
    else:
        if stage in NEEDS_MAP:
            env['ASSEMBLER_MAPPING_MANIFEST'] = str(mapping_result(stage, args, root, env, sources))
        if stage == 'chimera':
            taxonomy = env['GTDB_GENE_TAX_TSV']
            cmd = ['bash', str(ROOT / '07_chimera.sh'), str(taxonomy), str(root / '07_chimera')]
        elif stage == 'false-closure':
            for name in ('clipping_table', 'clipping_all_table', 'assembler_circularity'):
                need_result(name, args, root, env, stage)
            cmd = ['bash', str(ROOT / '04_false_closure.sh'), str(main_input), str(root / '04_false_closure')]
        elif stage == 'variants':
            cmd = ['bash', str(ROOT / '09_assembly_variants.sh'), str(main_input), str(root / '09_assembly_variants')]
        elif stage == 'phantom':
            cmd = ['bash', str(ROOT / '05_phantom.sh'), str(main_input), str(root / '05_phantom')]
        elif stage == 'repeats':
            cmd = ['bash', str(ROOT / '06_repeats.sh'), str(main_input), str(root / '06_repeats')]
        elif stage == 'false-duplication':
            need_result('clipping_table', args, root, env, stage)
            cmd = ['bash', str(ROOT / '08_false_duplication.sh'), str(main_input), str(root / '08_false_duplication')]
        elif stage == 'summary':
            if args.mapping_manifest:
                env['ASSEMBLER_MAPPING_MANIFEST'] = str(check_manifest(args.mapping_manifest))
            enabled = {name for name, var in FLAGS.items() if env.get(var, 'Y') == 'Y'}
            tables = {
                'clipping': ('clipping_table',),
                'chimera': ('chimera_table',),
                'circularity': ('circularity_table', 'circularity_audit'),
                'false-closure': ('false_closure_table', 'false_closure_audit'),
                'variants': ('variants_table',),
                'phantom': ('phantom_table',),
                'repeats': ('repeats_table', 'repeats_mag_table'),
                'false-duplication': ('duplication_table', 'duplication_clipping_table', 'duplication_clipping_mag_table'),
            }
            missing = []
            for name in sorted(enabled, key=list(FLAGS).index):
                for table in tables[name]:
                    if need_result(table, args, root, env, stage, optional=True) is None:
                        missing.append((name, FILE_OPTIONS[table][1], root / RESULTS[table]))
            if missing:
                done = [name for name in FLAGS if name in enabled and name not in {m[0] for m in missing}]
                first = {}
                for name, flag, path in missing:
                    first.setdefault(name, (flag, path))
                lines = [f'  {name}: {path} (or {flag} FILE)' for name, (flag, path) in first.items()]
                modules = ', '.join(first)
                tip = (f'\nOr summarize only the finished modules: MAGverify summary ... --include {",".join(done)}'
                       if done else '')
                raise UsageError('summary needs results that are not in this -o DIR yet:\n' + '\n'.join(lines)
                                 + f'\nRun first: MAGverify {modules.replace(", ", ", MAGverify ")}' + tip)
            cmd = ['bash', str(ROOT / '10_summary.sh'), str(main_input), str(root / '.cache' / 'summary')]
        elif stage == 'circularity':
            cmd = ['bash', str(ROOT / '02_circularity.sh'), env['MAG_DIR'], str(root / '02_circularity')]
        else:
            cmd = ['bash', str(ROOT / '03_clipping.sh'), str(main_input), str(root / '03_clipping')]
    print('Output: ' + str(root if stage == 'summary' else root / STAGES[stage][1]), flush=True)
    run(cmd, env, args.dry_run)
    if not args.dry_run:
        public_tables.publish(root)
        report_final_tables(stage, root)


def add_stage_parser(subparsers, stage):
    description = STAGES[stage][2]
    members = {stage}
    p = subparsers.add_parser(stage, description=description, help=description, add_help=False,
                              formatter_class=help_formatter)
    p.add_argument('-h', action='help', help=argparse.SUPPRESS)
    def handles(*names):
        return bool(members.intersection(names))
    if handles('map', 'clipping', 'chimera', 'circularity', 'false-closure', 'variants', 'phantom', 'false-duplication'):
        p.add_argument('-r', dest='reads', metavar='DIR',
                       help='Directory of HiFi reads (*.fastq.gz; *.fq.gz)')
    p.add_argument('-a', dest='source_assemblies', metavar='DIR',
                   help='Directory of assembled contigs (*.fa)')
    mag_stage = handles('clipping','chimera','circularity','false-closure','variants','phantom','repeats','false-duplication','summary')
    if mag_stage:
        p.add_argument('-g', dest='mag_dir', metavar='DIR',
                       help='Directory of MAGs (*.fa)')
    p.add_argument('-t', dest='threads', type=int, metavar='N', help='Number of threads')
    p.add_argument('-o', dest='output_directory', required=True, metavar='DIR', help='Directory of base output')
    if mag_stage and stage != 'circularity':
        p.add_argument('--mag-map', metavar='TSV', help=argparse.SUPPRESS)
    if handles('clipping','chimera','false-closure','variants','phantom','repeats','false-duplication','summary'):
        p.add_argument('--mapping-manifest', metavar='TSV', help=argparse.SUPPRESS)
    if handles('chimera','summary'):
        p.add_argument('--taxonomy', metavar='TSV', help='Gene taxonomy TSV (chimera output)')
    if handles('chimera','false-closure','summary'):
        p.add_argument('--scg-table', metavar='TSV', help='SCG TSV (prep-scg output)')
    if handles('chimera'):
        p.add_argument('--mag-taxonomy', metavar='TSV', help='MAG taxonomy TSV (prep-taxonomy output)')
    if handles('chimera','false-closure','false-duplication','summary'):
        p.add_argument('--clipping-table', metavar='TSV', help='Clipping site TSV (clipping output)')
    if handles('false-closure','summary'):
        p.add_argument('--clipping-all-table', metavar='TSV', help='All-region clipping TSV (clipping output)')
    if handles('false-closure','circularity','summary'):
        p.add_argument('--assembler-circularity', metavar='TSV',
                       help='Assembler circularity TSV (prep-circularity output)')
    if handles('summary'):
        for flag, message in [
            ('--circularity-table', 'Circularity MAG TSV (circularity output)'),
            ('--circularity-audit', 'Circularity audit TSV (circularity output)'),
            ('--false-closure-table', 'False closure MAG TSV (false-closure output)'),
            ('--false-closure-audit', 'False closure audit TSV (false-closure output)'),
            ('--chimera-table', 'Chimera contig TSV (chimera output)'),
            ('--variants-table', 'Variants contig TSV (variants output)'),
            ('--phantom-table', 'Phantom region TSV (phantom output)'),
            ('--repeats-table', 'Repeats contig TSV (repeats output)'),
            ('--repeats-mag-table', 'Repeats MAG TSV (repeats output)'),
            ('--duplication-table', 'False duplication contig TSV (false-duplication output)'),
            ('--duplication-clipping-table', 'False duplication clipping contig TSV (false-duplication output)'),
            ('--duplication-clipping-mag-table', 'False duplication clipping MAG TSV (false-duplication output)')]:
            p.add_argument(flag, metavar='TSV', help=message)
        p.add_argument('--include', metavar='MODULES', help='Completed modules (comma-separated)')
    if stage == 'chimera':
        p.add_argument('--database', metavar='PREFIX', help='MMseqs2 database prefix (prep-chimera output)')
    if handles('map', 'circularity', 'repeats', 'variants'):
        p.add_argument('-m', dest='memory_gb', type=int, metavar='N', help='Memory limit (GB)')
    if handles('variants'):
        p.add_argument('--memory-reserve-gb', type=int, metavar='N', help='Reserved free memory (GB, default: 64)')
    return p


DOCTOR_PROGRAMS = (
    ('bash', None), ('python3', None), ('awk', None), ('sed', None), ('xargs', None),
    ('sha256sum', None), ('minimap2', ['--version']), ('samtools', ['--version']),
    ('bedtools', ['--version']), ('blastn', ['-version']), ('makeblastdb', None),
    ('mmseqs', ['version']), ('pyrodigal', ['--version']), ('hmmsearch', ['-h']),
    ('meryl', ['--version']), ('anvi-profile', ['--version']),
    ('anvi-gen-contigs-database', None), ('anvi-run-hmms', None),
    ('anvi-script-find-misassemblies', None), ('anvi-script-reformat-fasta', None),
    ('anvi-script-reformat-bam', None), ('anvi-export-table', None),
)


def program_version(path, args):
    try:
        result = subprocess.run([path, *args], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return ''
    for line in (result.stdout + '\n' + result.stderr).splitlines():
        line = line.strip().lstrip('#').strip()
        if re.search(r'\d+\.\d+|\(v\d+', line):
            return re.sub(r'\s*\.{3,}\s*:?\s*', ': ', line)[:60]
    return ''


def doctor():
    """Check the installation: code location, required programs and their versions."""
    try:
        from importlib.metadata import version
        installed = version('MAGverify')
    except Exception:
        installed = 'not installed (running from the code folder)'
    checkout = ROOT.parent
    complete = all((checkout / name).is_file() for name in
                   ('magverify_cli.py', 'modules/lib_config.sh', 'scripts/core_magverify_cli.py'))
    print(f'{"MAGverify version":32s} {installed}')
    print(f'{"MAGverify path":32s} {checkout}' + ('' if complete else '  (INCOMPLETE: files missing)'))
    print(f'{"Python":32s} {sys.version.split()[0]}')
    print()
    print(f'{"Program":32s} {"Status":8s} Version')
    missing = []
    for name, version_args in DOCTOR_PROGRAMS:
        path = shutil.which(name)
        if not path:
            missing.append(name)
            print(f'{name:32s} {"MISSING":8s}')
            continue
        found = program_version(path, version_args) if version_args else ''
        print(f'{name:32s} {"OK":8s} {found}')
    print()
    if missing or not complete:
        if missing:
            print('Missing programs: ' + ', '.join(missing))
        print('Activate the MAGverify environment, or see the Installation section of README.md')
        return 1
    print(f'All {len(DOCTOR_PROGRAMS)} required programs found.')
    return 0


PIPELINE_STEPS = ('prep-circularity', 'prep-taxonomy', 'prep-scg', 'map', 'circularity', 'clipping',
                  'false-closure', 'phantom', 'repeats', 'false-duplication', 'variants', 'chimera', 'summary')


def add_run_parser(subparsers):
    p = subparsers.add_parser('end-to-end', formatter_class=help_formatter,
        help='Run all preparation and analysis modules in order',
        description='Run all preparation and analysis modules in order with one set of inputs. '
                    'Finished steps are reused, so an interrupted run can be restarted with the same command.',
        add_help=False)
    p.add_argument('-h', action='help', help=argparse.SUPPRESS)
    p.add_argument('-a', dest='source_assemblies', required=True, metavar='DIR',
                   help='Directory of assembled contigs (*.fa)')
    p.add_argument('-r', dest='reads', required=True, metavar='DIR',
                   help='Directory of HiFi reads (*.fastq.gz; *.fq.gz)')
    p.add_argument('-g', dest='mag_dir', required=True, metavar='DIR', help='Directory of MAGs (*.fa)')
    p.add_argument('-o', dest='output_directory', required=True, metavar='DIR', help='Directory of base output')
    p.add_argument('-t', dest='threads', type=int, metavar='N', help='Number of threads')
    p.add_argument('-m', dest='memory_gb', type=int, metavar='N', help='Memory limit (GB)')
    p.add_argument('-c', dest='sources', action='append', metavar='FILE',
                   help='Assembler circularity file for prep-circularity; repeat for each assembler '
                        '(not needed when 00_prep/assembler_circularity.MAG.tsv exists)')
    p.add_argument('-d', dest='gtdbtk_dir', metavar='DIR',
                   help='Directory of GTDB-Tk classify_wf results for prep-taxonomy '
                        '(not needed when 00_prep/gtdbtk_taxonomy.MAG.tsv exists)')
    p.add_argument('--database', metavar='PREFIX',
                   help='MMseqs2 database prefix for chimera (prep-chimera output; default: 00_prep/GTDB_R232)')
    p.add_argument('--skip', metavar='MODULES', default='',
                   help='Modules to leave out (comma-separated), e.g. --skip chimera,variants')
    p.add_argument('--memory-reserve-gb', type=int, metavar='N', help=argparse.SUPPRESS)
    return p


def pipeline_plan(args):
    """Return the (step, argv) list for MAGverify end-to-end, or raise UsageError before anything starts."""
    skip = {name.strip() for name in args.skip.split(',') if name.strip()}
    unknown = sorted(skip - set(PIPELINE_STEPS))
    if unknown:
        raise UsageError(f'--skip: unknown module(s) {", ".join(unknown)}; choose from {", ".join(PIPELINE_STEPS)}')
    if 'map' in skip:
        raise UsageError('--skip map is not possible: every analysis module reads the map result')
    dir_input(args.source_assemblies, '-a')
    dir_input(args.reads, '-r')
    dir_input(args.mag_dir, '-g')
    for path in args.sources or ():
        file_input(path, '-c')
    if args.gtdbtk_dir:
        dir_input(args.gtdbtk_dir, '-d')
    root = Path(args.output_directory).expanduser().resolve()
    if root.name in {'00_prep', '10_summary'} | {folder for _script, folder, _text in STAGES.values()}:
        root = root.parent
    have = lambda key: (root / RESULTS[key]).is_file() and (root / RESULTS[key]).stat().st_size > 0
    common = ['-o', str(root)] + (['-t', str(args.threads)] if args.threads else [])
    memory = ['-m', str(args.memory_gb)] if args.memory_gb else []
    if args.sources:
        pass
    elif have('assembler_circularity'):
        skip.add('prep-circularity')
    elif 'false-closure' not in skip:
        raise UsageError('false-closure needs the assembler circularity table: pass -c FILE for each assembler '
                         '(hifiasm-meta p_ctg.gfa, metaMDBG contigs.fasta.gz, metaFlye assembly_info.txt), '
                         'or --skip false-closure')
    else:
        skip.add('prep-circularity')
    if not args.gtdbtk_dir:
        skip.add('prep-taxonomy')
    if {'false-closure', 'chimera'} <= skip:
        skip.add('prep-scg')
    if 'chimera' not in skip and not have('taxonomy'):
        prefix = (Path(args.database).expanduser().resolve() if args.database else root / RESULTS['database'])
        if not Path(str(prefix) + '.dbtype').is_file():
            raise UsageError(f'chimera needs the MMseqs2 database ({prefix}): pass --database PREFIX '
                             '(MAGverify prep-chimera output), or --skip chimera')
    plan = []
    for step in PIPELINE_STEPS:
        if step in skip:
            continue
        if step == 'prep-circularity':
            argv = [step, '-g', args.mag_dir, '-o', str(root)] + [x for f in args.sources for x in ('-c', f)]
        elif step == 'prep-taxonomy':
            argv = [step, '-d', args.gtdbtk_dir, '-g', args.mag_dir, '-o', str(root)]
        elif step == 'prep-scg':
            argv = [step, '-a', args.source_assemblies] + common
        elif step == 'map':
            argv = [step, '-a', args.source_assemblies, '-r', args.reads] + common + memory
        elif step == 'summary':
            done = [name for name in FLAGS if name not in skip]
            argv = [step, '-g', args.mag_dir, '-o', str(root)]
            if len(done) < len(FLAGS):
                argv += ['--include', ','.join(done)]
        else:
            argv = [step, '-g', args.mag_dir] + common
            if step in ('circularity', 'repeats', 'variants'):
                argv += memory
            if step == 'variants' and args.memory_reserve_gb:
                argv += ['--memory-reserve-gb', str(args.memory_reserve_gb)]
            if step == 'chimera' and args.database:
                argv += ['--database', args.database]
        plan.append((step, argv))
    return plan


def pipeline_run(parser, args):
    plan = pipeline_plan(args)
    print('Steps: ' + ', '.join(step for step, _ in plan), flush=True)
    if args.dry_run:
        for step, argv in plan:
            print('MAGverify ' + shlex.join(argv), flush=True)
        return
    started = time.monotonic()
    for number, (step, argv) in enumerate(plan, 1):
        print(f'\n=== [{number}/{len(plan)}] MAGverify {step} ===', flush=True)
        dispatch(parser.parse_args(argv))
    print(f'\nAll steps finished. Total elapsed: {format_duration(time.monotonic() - started)}', flush=True)


def build_parser():
    parser = argparse.ArgumentParser(prog='MAGverify', description='PacBio HiFi MAG assembly quality control',
                                     formatter_class=help_formatter,
                                     add_help=False)
    parser.add_argument('-h', action='help', help=argparse.SUPPRESS)
    sub = parser.add_subparsers(dest='module', required=True)
    add_run_parser(sub)
    add_prep_parser(sub)
    for stage in STAGES:
        add_stage_parser(sub, stage)
    sub.add_parser('list', help='List available modules', add_help=False)
    sub.add_parser('doctor', help='Check the installation and required programs', add_help=False)
    return parser


def dispatch(args):
    args.dry_run = DRY_RUN
    if args.module == 'prep-chimera':
        prep_chimera_run(args)
    elif args.module == 'prep-taxonomy':
        prep_taxonomy_run(args)
    elif args.module == 'prep-circularity':
        prep_circularity_run(args)
    elif args.module == 'prep-scg':
        prep_scg_run(args)
    else:
        stage_run(args)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.dry_run = DRY_RUN
    if args.module == 'list':
        print('All modules in one command:')
        print('  end-to-end         Run the preparation and analysis modules below in order')
        print('Preparation:')
        print('  prep-chimera       Build MMseqs2 taxonomy database from local GTDB R232 files (once)')
        print('  prep-taxonomy      Build the MAG taxonomy table from GTDB-Tk results')
        print('  prep-circularity   Build the assembler circularity table')
        print("  prep-scg           Annotate single-copy genes with anvi'o")
        print('Analysis modules:')
        for name, (_, _, desc) in STAGES.items():
            print(f'  {name:18s} {desc}')
        return 0
    if args.module == 'doctor':
        return doctor()
    try:
        if args.module == 'end-to-end':
            pipeline_run(parser, args)
        else:
            dispatch(args)
    except UsageError as exc:
        parser.exit(2, f'ERROR: {exc}\nSee: MAGverify {args.module} -h\n')
    except subprocess.CalledProcessError as exc:
        parser.exit(exc.returncode or 1, f'ERROR: command failed ({exc.returncode}): {shlex.join(map(str, exc.cmd))}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
