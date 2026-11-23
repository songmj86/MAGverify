"""Regression tests for memory-efficient chimera depth and indexed HiFi windows."""
import io
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'scripts'))
import chimera_call as chimera  # noqa: E402


class ChimeraOptimizationTests(unittest.TestCase):
    def test_interval_index_matches_exhaustive_scan(self):
        rng = random.Random(723)
        for count in (1, 3, 100, 2000):
            windows = []
            for index in range(count):
                start = rng.randrange(0, 50000)
                end = start + rng.randrange(1, 2000)
                windows.append(chimera.ValidationWindow(index, 'contig', start, end, 0, 0, 0, 0, 0, 0))
            windows.sort(key=lambda x: (x.start0, x.index))
            starts = [x.start0 for x in windows]
            max_ends = []
            greatest = 0
            for item in windows:
                greatest = max(greatest, item.end0)
                max_ends.append(greatest)
            for _ in range(150):
                start = rng.randrange(0, 50000)
                end = start + rng.randrange(1, 3000)
                expected = {w.index for w in windows if w.start0 < end and w.end0 > start}
                observed = {w.index for w in chimera.overlapping_validation_windows(windows, starts, max_ends, start, end)}
                self.assertEqual(observed, expected)

    def test_validate_mode_matches_unindexed_evidence(self):
        anchor = 5
        lengths = {'ctg': 4000, 'other': 2000}
        positions = [(100, 105), (100, 105), (105, 120), (900, 901), (1090, 1100), (3000, 3000)]
        events = [dict(contig='ctg', transition_interval_start_1based=str(left),
                       transition_interval_end_1based=str(right)) for left, right in positions]
        events.append(dict(contig='other', transition_interval_start_1based='500',
                           transition_interval_end_1based='500'))
        reads = [
            ('r1', 'ctg', 50, '200M'),
            ('r2', 'ctg', 94, '10S30M'),
            ('r3', 'ctg', 90, '5M10D80M'),
            ('r4', 'ctg', 880, '100M'),
            ('r5', 'ctg', 1060, '50M'),
            ('r6', 'ctg', 2900, '200M'),
            ('r7', 'other', 480, '50M'),
            ('r8', 'ctg', 201, '100M'),
        ]
        sams = []
        for qname, contig, pos, cigar in reads:
            sams.append('\t'.join([qname, '0', contig, str(pos), '60', cigar, '*', '0', '0', '*', '*']))

        class FakeProc:
            def __init__(self):
                self.stdout = io.StringIO('\n'.join(sams) + '\n')
                self.stderr = io.StringIO('')
            def wait(self):
                return 0

        with patch.object(chimera.subprocess, 'Popen', return_value=FakeProc()) as fake:
            observed = chimera.validate_mode('fake.bam', events, list(range(len(events))), lengths, anchor, 2)
        self.assertIn('-L', fake.call_args.args[0])

        expected = {}
        for idx, event in enumerate(events):
            contig = event['contig']
            left, right, lo, hi = chimera.event_bounds(event, lengths[contig])
            start = max(0, left - anchor - 1)
            end = min(lengths[contig], right + anchor)
            expected[idx] = {'left': set(), 'right': set(), 'spanning': set(), 'clipped': set()}
            for qname, rcontig, pos1, cigar in reads:
                if rcontig != contig:
                    continue
                ref_start = pos1 - 1
                ops = chimera.cigar_ops(cigar)
                ref_end = ref_start + sum(n for n, op in ops if op in chimera.REF_OPS)
                if ref_end <= start or ref_start >= end:
                    continue
                runs = chimera.matched_runs(ref_start, ops)
                left_bp = chimera.matched_bp(runs, max(0, left - anchor), left)
                right_bp = chimera.matched_bp(runs, right - 1, min(lengths[contig], right - 1 + anchor))
                if left_bp >= anchor:
                    expected[idx]['left'].add(qname)
                if right_bp >= anchor:
                    expected[idx]['right'].add(qname)
                if left_bp >= anchor and right_bp >= anchor:
                    expected[idx]['spanning'].add(qname)
                if (chimera.terminal_clip(ops, True) > 0 and lo - 1 <= ref_start + 1 <= hi + 1) or (chimera.terminal_clip(ops, False) > 0 and lo - 1 <= ref_end <= hi + 1):
                    expected[idx]['clipped'].add(qname)
        for idx in expected:
            self.assertEqual(observed[idx], expected[idx], idx)

    def test_bedcov_streams_and_reports_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / 'samtools'
            script.write_text('''#!/usr/bin/env python3
import os, sys
if os.environ.get('MOCK_BEDCOV_FAIL'):
    print('failed bedcov', file=sys.stderr)
    sys.exit(7)
with open(sys.argv[2]) as bed:
    for line in bed:
        fields = line.rstrip('\\n').split('\\t')
        print('\\t'.join(fields + [str((int(fields[2])-int(fields[1])) * 12)]))
''')
            script.chmod(0o755)
            env = dict(os.environ, PATH=directory + os.pathsep + os.environ.get('PATH', ''))
            routes = {'c1': ('A', '/fakeA.bam'), 'c2': ('B', '/fakeB.bam')}
            intervals = [('c1', 10, 20, 'a'), ('c2', 100, 200, 'b'), ('c1', 22, 30, 'c')]
            with patch.dict(os.environ, env, clear=True):
                means = chimera.bedcov_means(routes, intervals, parallel_bam_jobs=2)
                self.assertEqual(means, {'a': 12.0, 'b': 12.0, 'c': 12.0})
                with patch.dict(os.environ, {'MOCK_BEDCOV_FAIL': '1'}):
                    with self.assertRaises(subprocess.CalledProcessError) as caught:
                        chimera.bedcov_means(routes, intervals, parallel_bam_jobs=1)
                    self.assertIn('failed bedcov', caught.exception.stderr)


if __name__ == '__main__':
    unittest.main()
