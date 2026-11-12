#!/usr/bin/env python3
"""Command-line entry point for MAGverify."""
import argparse
import csv
import json
import os
import re
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent.parent / "modules"
SCRIPTS_DIR = Path(__file__).resolve().parent
STAGES = {
    'map': ('01_mapping.sh', '01_mapping', 'Map HiFi reads to source assemblies'),
    'circularity': ('02_circularity.sh', '02_circularity', 'Validate circular genome topology'),
    'clipping': ('03_clipping.sh', '03_clipping', 'Assess read-clipping sites'),
    'false-closure': ('04_false_closure.sh', '04_false_closure', 'Assess false genome closure'),
    'phantom': ('05_phantom.sh', '05_phantom', 'Detect phantom sequences'),
    'repeats': ('06_repeats.sh', '06_repeats', 'Characterize repeat regions'),
    'chimera': ('07_chimera.sh', '07_chimera', 'Assess gene-level taxonomic chimeras'),
    'false-duplication': ('08_false_duplication.sh', '08_false_duplication', 'Assess false duplication'),
    'variants': ('09_assembly_variants.sh', '09_assembly_variants', 'Assess assembly variants'),
    'summary': ('10_summary.sh', '10_summary', 'Summarize and validate QC results'),
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
    'taxonomy': '00_prep/gene_taxonomy/GTDB_R232/gene_taxonomy.mmseqs2.GTDB_R232.all_contigs.tsv',
    'mag_map': '00_prep/MAG_contig_map.tsv',
    'mapping_manifest': '01_mapping/assembler_mapping_outputs.tsv',
    'clipping_table': '03_clipping/clipping.filtered.tsv',
    'clipping_all_table': '03_clipping/clipping.all_regions.filtered.tsv',
    'scg_table': '00_prep/SCG/SCG_from_anvio.original_contig.tsv',
    'circularity_table': '02_circularity/MAG_circularity.final.tsv',
    'circularity_audit': '02_circularity/MAG_circularity.audit.tsv',
    'false_closure_table': '04_false_closure/MAG_false_closure.final.tsv',
    'false_closure_audit': '04_false_closure/MAG_false_closure.audit.tsv',
    'chimera_table': '07_chimera/GTDB_R232/clipping_defined/contig_taxonomic_chimera.GTDB_R232.final.tsv',
    'variants_table': '09_assembly_variants/contig_assembly_allele_support.final.tsv',
    'phantom_table': '05_phantom/phantom_sequences.final.tsv',
    'repeats_table': '06_repeats/contig_repeat_breadth.final.tsv',
    'repeats_mag_table': '06_repeats/MAG_repeat_breadth.final.tsv',
    'duplication_table': '08_false_duplication/false_duplication.by_contig.tsv',
    'duplication_clipping_table': '08_false_duplication/validation/clipping/false_duplication.by_contig.clipping.tsv',
    'duplication_clipping_mag_table': '08_false_duplication/validation/clipping/false_duplication.by_MAG.clipping.tsv',
}
# Stable, short filenames for the primary user-facing tables. Internal filenames
# remain unchanged to preserve inputs and provenance across modules.
PUBLIC_RESULTS = {
    'map': {
        'mapping_manifest.tsv': '01_mapping/assembler_mapping_outputs.tsv',
    },
    'clipping': {
        'clipping_sites.tsv': '03_clipping/clipping.filtered.tsv',
        'clipping_mag.tsv': '03_clipping/MAG_clipping.final.tsv',
    },
    'chimera': {
        'gene_taxonomy.tsv': '00_prep/gene_taxonomy/GTDB_R232/gene_taxonomy.mmseqs2.GTDB_R232.all_contigs.tsv',
        'scg.tsv': '00_prep/SCG/SCG_from_anvio.original_contig.tsv',
        'chimera_contig.tsv': '07_chimera/GTDB_R232/clipping_defined/contig_taxonomic_chimera.GTDB_R232.final.tsv',
        'chimera_mag.tsv': '07_chimera/GTDB_R232/clipping_defined/MAG_taxonomic_chimera.GTDB_R232.final.tsv',
    },
    'circularity': {
        'circularity_mag.tsv': '02_circularity/MAG_circularity.final.tsv',
    },
    'false-closure': {
        'false_closure_mag.tsv': '04_false_closure/MAG_false_closure.final.tsv',
        'false_closure_contig.tsv': '04_false_closure/contig_false_closure.final.tsv',
    },
    'variants': {
        'variants_contig.tsv': '09_assembly_variants/contig_assembly_allele_support.final.tsv',
        'variants_mag.tsv': '09_assembly_variants/assembly_allele_support.by_MAG.tsv',
    },
    'phantom': {
        'phantom_contig.tsv': '05_phantom/phantom_sequences.final.tsv',
        'phantom_mag.tsv': '05_phantom/MAG_phantom_summary.final.tsv',
    },
    'repeats': {
        'repeats_contig.tsv': '06_repeats/contig_repeat_breadth.final.tsv',
        'repeats_mag.tsv': '06_repeats/MAG_repeat_breadth.final.tsv',
    },
    'false-duplication': {
        'false_duplication_contig.tsv': '08_false_duplication/false_duplication.by_contig.tsv',
        'false_duplication_mag.tsv': '08_false_duplication/false_duplication.by_MAG.tsv',
    },
    'summary': {
        'summary_mag.tsv': '10_summary/integrated_QC_summary.by_MAG.tsv',
        'summary_contig.tsv': '10_summary/integrated_QC_summary.MAG_contigs.tsv',
    },
}


def publish_results(stage, root):
    """Expose final tables with simple names without duplicating large files."""
    names = PUBLIC_RESULTS.get(stage, {})
    available = [(name, root / relative) for name, relative in names.items()
                 if (root / relative).is_file() and (root / relative).stat().st_size > 0]
    if not available:
        return
    public_dir = root / 'results'
    public_dir.mkdir(parents=True, exist_ok=True)
    for name, original in available:
        short_path = public_dir / name
        relative_target = os.path.relpath(original, public_dir)
        if short_path.is_symlink():
            if short_path.resolve() == original.resolve():
                continue
            short_path.unlink()
        elif short_path.exists():
            raise UsageError(f'Output name already exists and is not a managed link: {short_path}')
        short_path.symlink_to(relative_target)
    print('Final tables:', flush=True)
    for name, _ in available:
        print('  ' + str(public_dir / name), flush=True)


FLAGS = {
    'clipping': 'RUN_CLIPPING', 'chimera': 'RUN_TAXONOMIC_CHIMERA',
    'circularity': 'RUN_CIRCULARITY', 'false-closure': 'RUN_FALSE_CLOSURE',
    'variants': 'RUN_ASSEMBLY_VARIANTS', 'phantom': 'RUN_PHANTOM',
    'repeats': 'RUN_REPEAT', 'false-duplication': 'RUN_FALSE_DUPLICATION',
}
NEEDS_MAP = {'clipping', 'chimera', 'false-closure', 'variants', 'phantom', 'false-duplication'}
NEEDS_MAG_MAP = {'chimera', 'false-closure', 'variants', 'phantom', 'repeats', 'false-duplication', 'summary'}
NEEDS_SOURCES = {'map', 'circularity', 'repeats'}


class UsageError(Exception):
    pass


def file_input(value, name):
    if not value:
        raise UsageError(f'Missing {name}. Run MAGverify MODULE --help for input options.')
    path = Path(value).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise UsageError(f'{name} must be an existing non-empty file: {path}')
    return path


def dir_input(value, name):
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise UsageError(f'{name} must be an existing directory: {path}')
    return path


def assembly_fasta(value):
    path = file_input(value, '--input')
    if path.suffix.lower() not in ('.fa', '.fasta', '.fna'):
        raise UsageError('--input must be an uncompressed assembly FASTA (.fa, .fna, .fasta), not a study directory')
    return path


def fasta_paths(directory, label):
    parent = dir_input(directory, label)
    paths = sorted(p for p in parent.iterdir() if p.is_file() and p.suffix.lower() in ('.fa', '.fna', '.fasta'))
    if not paths:
        raise UsageError(f'No FASTA files found in {parent}; check the directory or provide individual FASTAs')
    return paths


def validate_fasta_paths(paths, option):
    resolved = [file_input(p, option) for p in paths]
    if any(p.suffix.lower() not in ('.fa', '.fna', '.fasta') for p in resolved):
        raise UsageError(f'{option} accepts .fa, .fasta and .fna files')
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
    if getattr(args, 'source_assembly', None) and getattr(args, 'source_assemblies', None):
        raise UsageError('Choose repeated --source-assembly FILE options OR --assemblies-dir DIR, not both')
    if getattr(args, 'source_assembly', None):
        specified = []
        for value in args.source_assembly:
            label = None
            raw = value
            if '=' in value:
                label, raw = value.split('=', 1)
                if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', label):
                    raise UsageError('--source-assembly labels must use letters, digits, dots, dashes or underscores')
            source = file_input(raw, '--source-assembly')
            if source.suffix.lower() not in ('.fa', '.fasta', '.fna'):
                raise UsageError('--source-assembly accepts .fa, .fna or .fasta')
            specified.append((label, source))
        if len({path for _, path in specified}) != len(specified):
            raise UsageError('The same source FASTA was supplied more than once')
        names = []
        for label, source in specified:
            name = f'{label}{source.suffix}' if label else source.name
            if not re.fullmatch(r'[A-Za-z0-9_.-]+\.(fa|fasta|fna)', name, re.IGNORECASE):
                name = re.sub(r'[^A-Za-z0-9_.-]+', '_', source.stem) + source.suffix
            if name in names and not label:
                name = re.sub(r'[^A-Za-z0-9_.-]+', '_', source.parent.name) + '_' + name
            if name in names:
                raise UsageError('Source assembly names are ambiguous. Use --source-assembly NAME=/path/to/assembly.fasta')
            names.append(name)
        if len({Path(name).stem for name in names}) != len(names):
            raise UsageError('Source assembly names must have unique identifiers. Use NAME=/path/to/assembly.fasta')
        named = list(zip(names, [path for _, path in specified]))
        parent = link_inputs(named, output_root / '00_prep' / 'source_assembly_links', args.dry_run)
    elif getattr(args, 'source_assemblies', None):
        paths = validate_fasta_paths(fasta_paths(args.source_assemblies, '--assemblies-dir'), '--assemblies-dir')
        names = [p.name for p in paths]
        parent = dir_input(args.source_assemblies, '--assemblies-dir')
    else:
        return None
    print(f'Source assemblies ({len(names)}): ' + ', '.join(names), flush=True)
    return parent, ' '.join(names)


def mag_inputs(args, output_root):
    if getattr(args, 'mag_fasta', None) and getattr(args, 'mag_dir', None):
        raise UsageError('Choose repeated --mag-fasta FILE options OR --mag-dir DIR, not both')
    if getattr(args, 'mag_fasta', None):
        paths = validate_fasta_paths(args.mag_fasta, '--mag-fasta')
        parent = link_inputs([(path.name, path) for path in paths], output_root / '00_prep' / 'mag_fasta_links', args.dry_run)
    elif getattr(args, 'mag_dir', None):
        paths = validate_fasta_paths(fasta_paths(args.mag_dir, '--mag-dir'), '--mag-dir')
        parent = dir_input(args.mag_dir, '--mag-dir')
    else:
        return None
    print(f'MAG FASTAs ({len(paths)}): ' + ', '.join(path.name for path in paths), flush=True)
    return parent


def run(cmd, env, dry_run):
    print('$ ' + shlex.join(map(str, cmd)), flush=True)
    if not dry_run:
        subprocess.run(cmd, cwd=ROOT, env=env, check=True)


def check_programs(stage, mapping=False, env=None):
    core = {'bash', 'python3'}
    if stage in ('map', 'clipping'):
        core |= {'samtools', 'sha256sum'}
    if stage == 'map' or mapping:
        core.add('minimap2')
    if stage == 'clipping':
        core.add('anvi-script-find-misassemblies')
    if stage == 'prep-chimera':
        core.update({'wget', 'mmseqs', 'awk', 'sed'})
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
        raise UsageError('Missing programs: ' + ', '.join(missing) + '. See docs/INSTALL.md')


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
            file_input(str(bam) + '.inputs.sha256', f'BAM mapping provenance in {path}')
            count += 1
        if not count:
            raise UsageError(f'--mapping-manifest is empty: {path}')
    return path


def assemble_combined(paths, output, dry_run):
    """Build a FASTA spanning the original source assemblies, retaining contig IDs."""
    sources = [file_input(p, 'assembler FASTA') for p in paths]
    if not sources:
        raise UsageError('No assembler FASTAs available to build the combined assembly')
    fingerprint = [dict(path=str(p.resolve()), size=p.stat().st_size,
                        mtime_ns=p.stat().st_mtime_ns) for p in sources]
    metadata = output.with_suffix(output.suffix + '.sources.json')
    if output.is_file() and output.stat().st_size and metadata.is_file():
        try:
            if json.loads(metadata.read_text(encoding='utf-8')) == fingerprint:
                print(f'Combined assembly: {output} (reused)', flush=True)
                return output
        except (ValueError, OSError):
            pass
    print(f'Combined assembly: {output} ({len(sources)} source FASTA(s))', flush=True)
    if dry_run:
        return output
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix='.combined_', suffix='.fasta', dir=output.parent)
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
                                raise UsageError(f'Duplicate contig ID {identifier!r} across source assemblies. '
                                                 'IDs must be unique; inspect your assembler FASTAs.')
                            identifiers.add(identifier)
                            seen_header = True
                        elif line.strip() and not seen_header:
                            raise UsageError(f'Sequence data appears before a FASTA header in {source}')
                        dest.write(line if line.endswith(b'\n') else line + b'\n')
                if not seen_header:
                    raise UsageError(f'No FASTA records in {source}')
        os.replace(temp, output)
        metadata.write_text(json.dumps(fingerprint, indent=2) + '\n', encoding='utf-8')
        # A previously generated .fai cannot be reused after rebuilding the FASTA.
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
    directory = link_inputs(named, root / '00_prep' / 'source_assembly_links', args.dry_run)
    return directory, ' '.join(name for name, _ in named)


def defaults_for(stage, args):
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    threads = args.threads if args.threads is not None else min(100, cpus)
    if threads < 1:
        raise UsageError('--threads must be a positive integer')
    if threads > cpus:
        raise UsageError(f'--threads={threads} exceeds the available CPUs ({cpus})')
    root = Path(args.output_directory).expanduser().resolve()
    env = os.environ.copy()
    # No inherited project configuration or stale dataset-specific variables.
    for key in list(env):
        if key in {'OUT', 'READS', 'ASSEMBLY_FA', 'MAG_DIR', 'MAG_CONTIG_MAP',
                   'ASSEMBLER_ASSEMBLY_DIR', 'ASSEMBLER_FASTAS', 'ASSEMBLER_MAPPING_MANIFEST',
                   'MAPPING_DIR', 'CLIPPING_DIR', 'PREP_DIR', 'PER_ASSEMBLER_MANIFEST',
                   'PER_ASSEMBLER_MAP_ROOT', 'HIFI_MAPPING_DIR', 'SCG_TSV', 'SCG_OUT_DIR',
                   'CIRCULARITY_TSV', 'MAG_GTDB_TAXONOMY_TSV', 'GTDB_GENE_TAX_TSV',
                   'CLIP_FINAL', 'CLIP_ALL_FINAL', 'GTDB_CHIMERA_REPORT_CONTIG',
                   'CIRCULARITY_ONLY_FINAL', 'CIRCULARITY_ONLY_AUDIT', 'FALSE_CLOSURE_FINAL',
                   'FALSE_CLOSURE_AUDIT', 'ALLELE_CONTIG', 'PHANTOM_FINAL', 'REPEAT_SUM',
                   'REPEAT_MAG_SUM', 'FALSE_DUP_CONTIG_SUM', 'FALSE_DUP_CLIPPING_CONTIG_SUM',
                   'FALSE_DUP_CLIPPING_MAG_SUM'} or key.startswith((
                       'MAPPING_INPUT_', 'CLIPPING_INPUT_', 'CHIMERA_INPUT_', 'CIRCULARITY_INPUT_',
                       'FALSE_CLOSURE_INPUT_', 'VARIANT_INPUT_', 'PHANTOM_INPUT_', 'REPEAT_INPUT_',
                       'FALSE_DUP_INPUT_', 'SUMMARY_INPUT_')):
            env.pop(key, None)
    env.update(OUT=str(root), THREADS=str(threads), QC_SKIP_LOCAL_CONFIG='1')
    return root, env


def prepare_map(args, root, env, main_input):
    map_file = root / RESULTS['mag_map']
    if getattr(args, 'mag_map', None):
        return file_input(args.mag_map, '--mag-map')
    if 'MAG_DIR' not in env:
        if map_file.is_file():
            return file_input(map_file, 'MAG map in output directory')
        raise UsageError('MAG membership is needed: pass --mag-fasta FILE (repeatable), --mag-dir DIR, or --mag-map TSV')
    if not main_input:
        raise UsageError('Preparing a MAG map needs --input assembly FASTA')
    if args.dry_run:
        run(['samtools', 'faidx', str(main_input)], env, True)
        run(['python3', str(SCRIPTS_DIR / 'prep_build_mag_contig_map.py'), '--mag-dir', env['MAG_DIR'],
             '--extensions', 'fa fasta fna', '--assembly-fai', str(main_input) + '.fai',
             '--output', str(map_file)], env, True)
        return map_file
    if not shutil.which('samtools'):
        raise UsageError('samtools is necessary to index the assembly and prepare MAG membership')
    run(['samtools', 'faidx', str(main_input)], env, False)
    map_file.parent.mkdir(parents=True, exist_ok=True)
    run(['python3', str(SCRIPTS_DIR / 'prep_build_mag_contig_map.py'), '--mag-dir', env['MAG_DIR'],
         '--extensions', 'fa fasta fna', '--assembly-fai', str(main_input) + '.fai',
         '--output', str(map_file)], env, False)
    return file_input(map_file, 'generated MAG map')


def mapping_result(stage, args, root, env, sources):
    if args.mapping_manifest:
        return check_manifest(args.mapping_manifest)
    generated = root / RESULTS['mapping_manifest']
    if generated.is_file():
        manifest = check_manifest(generated)
        # An existing mapping belongs to a particular set of source FASTAs.
        # Refuse to silently reuse BAMs for different assembly inputs.
        if sources:
            requested = {p.resolve() for p in source_paths(sources)}
            recorded = {p.resolve() for p in source_fastas_from_manifest(manifest)}
            if requested != recorded:
                raise UsageError('Source assemblies differ from the existing mapping manifest. '
                                 'Rerun MAGverify map with the new inputs or use another --output-dir.')
        return manifest
    if sources and 'READS' not in env:
        raise UsageError('Run MAGverify map with --reads first, or supply --reads to map automatically')
    if 'READS' in env and sources:
        if not args.dry_run:
            check_programs('map', env=env)
        run(['bash', str(ROOT / '01_mapping.sh'), env['READS'], str(root / '01_mapping')], env, args.dry_run)
        return generated if args.dry_run else check_manifest(generated)
    if getattr(args, 'reads', None) and not sources and not generated.is_file():
        raise UsageError('Provide one or more --source-assembly FILE paths with --reads, or --mapping-manifest')
    raise UsageError(f'{stage} needs an existing map result in --output-dir, '
                     'or --reads FASTQ and --assemblies-dir DIR (or repeated --source-assembly FASTA).')


def choose_inputs(args, root, env):
    stage = args.module
    main_input = None
    if getattr(args, 'input', None):
        main_input = assembly_fasta(args.input)
        env['ASSEMBLY_FA'] = str(main_input)
    if getattr(args, 'reads', None):
        env['READS'] = str(file_input(args.reads, '--reads'))
    sources = source_inputs(args, root)
    if not sources and stage != 'map':
        sources = source_links_from_manifest(args, root)
    if sources:
        env['ASSEMBLER_ASSEMBLY_DIR'], env['ASSEMBLER_FASTAS'] = str(sources[0]), sources[1]
    if stage not in ('map', 'circularity') and not main_input:
        if sources:
            main_input = assemble_combined(source_paths(sources),
                                           root / '00_prep' / 'combined_assembly.fasta', args.dry_run)
        else:
            stored = root / '00_prep' / 'combined_assembly.fasta'
            if stored.is_file() and stored.stat().st_size:
                main_input = stored
                print(f'Combined assembly: {stored} (existing)', flush=True)
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
    # The active Conda environment supplies anvi'o and the other tools.
    # This also resolves bin/python paths used by the variants and clipping modules.
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
    if stage == 'variants' and args.memory_gb is not None and args.memory_gb < 1:
        raise UsageError('--memory-gb must be a positive integer')
    if stage == 'variants' and args.memory_gb:
        env['ASSEMBLY_VARIANT_MEMORY_BUDGET_GB'] = str(args.memory_gb)
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
    raise UsageError(f'{stage} needs {flag} FILE or an existing result at {fallback}')


def prep_chimera_run(args):
    """Prepare the GTDB R232 MMseqs2 database in OUTPUT/0_prep."""
    root, env = defaults_for('prep-chimera', args)
    script = ROOT / '00_prep_chimera.sh'
    if not script.is_file():
        raise UsageError(f'Missing installed GTDB preparation module: {script}')
    if not args.dry_run:
        check_programs('prep-chimera', env=env)
    prep = root / '0_prep'
    print('GTDB R232 database output: ' + str(prep / 'GTDB_R232'), flush=True)
    run(['bash', str(script), str(prep), env['THREADS']], env, args.dry_run)


def add_prep_parser(subparsers):
    """A setup command, separate from the ten analysis stages."""
    p = subparsers.add_parser('prep-chimera',
        help='Download GTDB R232 and build the MMseqs2 taxonomy database',
        description='Prepare the fixed-release GTDB R232 taxonomy database for the chimera module')
    p.add_argument('-o', '--output-dir', '--output-directory', dest='output_directory',
                   required=True, metavar='DIR',
                   help='Shared MAGverify output directory; database is created in DIR/0_prep')
    p.add_argument('-t', '--threads', type=int, metavar='N',
                   help='CPU threads (default up to 100, limited by detected CPUs)')
    p.add_argument('--dry-run', action='store_true',
                   help='Print preparation command without downloading or creating files')
    return p

def stage_run(args):
    stage = args.module
    root, env = defaults_for(stage, args)
    main_input, sources = choose_inputs(args, root, env)
    if stage not in ('map', 'circularity') and not main_input:
        raise UsageError(f'{stage} needs --assemblies-dir DIR (to build its assembly FASTA automatically), '
                         '--mapping-manifest TSV, or an existing --input FASTA')
    if stage in NEEDS_SOURCES and not sources:
        raise UsageError(f'{stage} needs --assemblies-dir DIR (original assembler FASTAs), '
                         'one or more --source-assembly FASTA options, or a mapping manifest')
    if stage == 'map' and not args.reads:
        raise UsageError('map needs --reads /path/to/HiFi.fastq.gz')
    if stage == 'circularity' and ('MAG_DIR' not in env or not args.reads):
        raise UsageError('circularity needs --reads plus --mag-fasta FILE(s) or --mag-dir DIR')
    if stage == 'clipping' and 'MAG_DIR' not in env:
        raise UsageError('clipping needs --mag-dir DIR (one FASTA per MAG); --mag-fasta is also supported')
    # Check stage-specific input tables before starting an expensive mapping run.
    if stage == 'chimera':
        # Chimera generates gene taxonomy itself when it is not already available.
        existing_taxonomy = need_result('taxonomy', args, root, env, stage, optional=True)
        if not existing_taxonomy:
            if not args.database:
                raise UsageError('chimera needs --database MMSEQS_PREFIX to generate gene taxonomy, '
                                 'or --taxonomy TSV to reuse a previous result')
            prefix = Path(args.database).expanduser().resolve()
            if not Path(str(prefix) + '.dbtype').is_file():
                raise UsageError(f'--database must be a MMseqs2 database prefix with .dbtype: {prefix}')
            env['MMSEQS_TAX_DB'] = str(prefix)
            env['GTDB_GENE_TAX_TSV'] = str(root / RESULTS['taxonomy'])
        need_result('clipping_table', args, root, env, stage)
    if stage == 'false-closure':
        for table in ('clipping_table', 'clipping_all_table', 'assembler_circularity'):
            need_result(table, args, root, env, stage)
    if stage == 'false-duplication':
        need_result('clipping_table', args, root, env, stage)
    if stage in ('variants', 'phantom') and not args.reads:
        raise UsageError(f'{stage} needs --reads FASTQ and the source mapping evidence')
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
            if not args.reads:
                raise UsageError('variants needs --reads FASTQ as well as a mapping manifest')
            cmd = ['bash', str(ROOT / '09_assembly_variants.sh'), str(main_input), str(root / '09_assembly_variants')]
        elif stage == 'phantom':
            if not args.reads:
                raise UsageError('phantom needs --reads FASTQ as well as a mapping manifest')
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
            for name in enabled:
                for table in tables[name]:
                    need_result(table, args, root, env, stage)
            cmd = ['bash', str(ROOT / '10_summary.sh'), str(main_input), str(root / '10_summary')]
        elif stage == 'circularity':
            cmd = ['bash', str(ROOT / '02_circularity.sh'), env['MAG_DIR'], str(root / '02_circularity')]
        else:
            cmd = ['bash', str(ROOT / '03_clipping.sh'), str(main_input), str(root / '03_clipping')]
    print('Output: ' + str(root / STAGES[stage][1]), flush=True)
    run(cmd, env, args.dry_run)
    if not args.dry_run:
        publish_results(stage, root)


def add_stage_parser(subparsers, stage):
    description = STAGES[stage][2]
    members = {stage}
    p = subparsers.add_parser(stage, description=description, help=description)
    def handles(*names):
        return bool(members.intersection(names))
    if stage not in ('map', 'circularity'):
        p.add_argument('-i', '--input', metavar='FASTA', help='Optional: existing combined assembly FASTA (otherwise generated automatically)')
    p.add_argument('-o', '--output-dir', '--output-directory', dest='output_directory', required=True, metavar='DIR', help='Output directory (reuse the same path across modules)')
    p.add_argument('-t', '--threads', type=int, metavar='N', help='CPU threads (100 recommended; default up to 100, limited by detected CPUs)')
    p.add_argument('--dry-run', action='store_true', help='Print the analysis command without running it')
    if handles('map', 'clipping', 'chimera', 'circularity', 'false-closure', 'variants', 'phantom', 'false-duplication'):
        p.add_argument('-r', '--reads', metavar='FASTQ', help='Original PacBio HiFi reads (FASTQ or FASTQ.GZ)')
    p.add_argument('--source-assembly', action='append', metavar='FASTA', help='Original assembler FASTA; repeat if needed. Accepts NAME=/path/to/file.fa')
    p.add_argument('-a', '--assemblies-dir', '--assemblies', '--source-assemblies', dest='source_assemblies', metavar='DIR', help='Directory of original assembler FASTA files (legacy aliases: --assemblies, --source-assemblies)')
    if handles('clipping','chimera','circularity','false-closure','variants','phantom','repeats','false-duplication','summary'):
        p.add_argument('--mag-fasta', action='append', metavar='FASTA', help='Alternative to --mag-dir: repeat for individual MAG FASTAs')
        p.add_argument('-m', '--mag-dir', metavar='DIR', help='Directory with one FASTA per MAG')
        p.add_argument('--mag-map', metavar='TSV', help='Existing MAG/contig table; otherwise generated')
    if handles('map','clipping','chimera','circularity','false-closure','variants','phantom','repeats','false-duplication','summary'):
        p.add_argument('--mapping-manifest', metavar='TSV', help='Reuse an existing mapping manifest')
    if handles('chimera','summary'):
        p.add_argument('--taxonomy', metavar='TSV', help='Gene taxonomy TSV')
    if handles('chimera','false-closure','summary'):
        p.add_argument('--scg-table', metavar='TSV', help='Single-copy marker gene table')
    if handles('chimera'):
        p.add_argument('--mag-taxonomy', metavar='TSV', help='MAG GTDB classification table')
    if handles('chimera','false-closure','false-duplication','summary'):
        p.add_argument('--clipping-table', metavar='TSV', help='Clipping classification TSV')
    if handles('false-closure','summary'):
        p.add_argument('--clipping-all-table', metavar='TSV', help='All-region clipping classification TSV')
    if handles('false-closure','circularity','summary'):
        p.add_argument('--assembler-circularity', metavar='TSV', help='Assembler-reported circularity TSV (necessary for false closure)')
    if handles('summary'):
        for flag, message in [
            ('--circularity-table','Circularity final TSV'),('--circularity-audit','Circularity audit TSV'),
            ('--false-closure-table','False closure TSV'),('--false-closure-audit','False closure audit TSV'),
            ('--chimera-table','Chimera final TSV'),('--variants-table','Assembly variants TSV'),
            ('--phantom-table','Phantom TSV'),('--repeats-table','Contig repeat TSV'),
            ('--repeats-mag-table','MAG repeat TSV'),('--duplication-table','False duplication TSV'),
            ('--duplication-clipping-table','False duplication clipping contig TSV'),
            ('--duplication-clipping-mag-table','False duplication clipping MAG TSV')]:
            p.add_argument(flag, metavar='TSV', help=message)
        p.add_argument('--include', metavar='MODULES', help='Comma-separated completed assessments for summary and validation')
    if stage == 'chimera':
        p.add_argument('--database', metavar='PREFIX', help='MMseqs2 gene taxonomy database prefix; omit if using --taxonomy TSV or a previous result')
    if handles('variants'):
        p.add_argument('--memory-gb', type=int, metavar='GB', help='Memory budget for assembly-variant profiling')
        p.add_argument('--memory-reserve-gb', type=int, metavar='GB', help='Reserve free memory (default 64 GiB)')
    return p


def main(argv=None):
    parser = argparse.ArgumentParser(prog='MAGverify', description='PacBio HiFi MAG assembly quality control')
    sub = parser.add_subparsers(dest='module', required=True)
    add_prep_parser(sub)
    for stage in STAGES:
        add_stage_parser(sub, stage)
    sub.add_parser('list', help='List available modules')
    sub.add_parser('doctor', help='List core bioinformatics commands available on PATH')
    args = parser.parse_args(argv)
    if args.module == 'list':
        print('Preparation (run once):')
        print('  prep-chimera       Download GTDB R232 and build MMseqs2 taxonomy database')
        print('Analysis modules (run separately):')
        for name, (_, _, desc) in STAGES.items():
            print(f'  {name:18s} {desc}')
        return 0
    if args.module == 'doctor':
        for cmd in ['bash', 'python3', 'minimap2', 'samtools', 'mmseqs', 'anvi-profile',
                    'anvi-script-find-misassemblies', 'blastn', 'makeblastdb', 'meryl', 'bedtools']:
            print(f'{cmd:34s} {shutil.which(cmd) or "not found"}')
        return 0
    try:
        if args.module == 'prep-chimera':
            prep_chimera_run(args)
        else:
            stage_run(args)
    except UsageError as exc:
        parser.exit(2, f'ERROR: {exc}\nSee: MAGverify {args.module} --help\n')
    except subprocess.CalledProcessError as exc:
        parser.exit(exc.returncode or 1, f'ERROR: command failed ({exc.returncode}): {shlex.join(map(str, exc.cmd))}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
