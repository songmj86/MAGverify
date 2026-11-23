"""Keep primary and supplementary alignments, omitting secondary in samtools."""
import io
from pathlib import Path
import sys
import unittest
from collections import defaultdict
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'scripts'))
import false_duplication as fd  # noqa: E402


class FilterTests(unittest.TestCase):
    def test_filter_keeps_primary_supplementary_only(self):
        rows = '\n'.join([
            '\t'.join([qname, str(flag), 'ctg', '1', '60', '100M', '*', '0', '0', '*', '*'])
            for qname, flag in [('primary', 0), ('supp', 2048)]
        ]) + '\n'

        class FakeProcess:
            def __init__(self):
                self.stdout = io.StringIO(rows)
            def wait(self):
                return 0

        evidence = {('ctg', 50): fd.BoundaryEvidence()}
        spanning = defaultdict(lambda: defaultdict(set))
        out = []
        with patch.object(fd.subprocess, 'Popen', return_value=FakeProcess()) as popen:
            fd.process_bam('test.bam', 'samtools', 5, {'ctg': [50]}, evidence, spanning, out)
        self.assertEqual(popen.call_args.args[0][2:4], ['-F', '260'])
        self.assertEqual(evidence['ctg', 50].spanning_reads, {'primary', 'supp'})
        self.assertEqual(spanning['ctg', 50]['primary'], {'primary'})
        self.assertEqual(spanning['ctg', 50]['supplementary'], {'supp'})


if __name__ == '__main__':
    unittest.main()
