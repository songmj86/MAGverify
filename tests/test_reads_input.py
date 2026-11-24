"""-r accepts a HiFi FASTQ file or a directory of them."""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from scripts.core_magverify_cli import UsageError, reads_input, reads_path_input  # noqa: E402

READ = '@r1\nACGT\n+\nIIII\n'


class ReadsInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)
        self.reads = self.dir / 'reads'
        self.reads.mkdir()
        (self.reads / 'a.fastq').write_text(READ)
        self.out = self.dir / 'out'

    def test_single_file_is_used_as_is(self):
        path = reads_input(self.reads / 'a.fastq', self.out, True)
        self.assertEqual(path, self.reads / 'a.fastq')
        self.assertFalse((self.out / '00_prep').exists())

    def test_directory_with_one_file_is_used_as_is(self):
        self.assertEqual(reads_input(self.reads, self.out, True), self.reads / 'a.fastq')

    def test_directory_with_several_files_is_joined(self):
        (self.reads / 'b.fastq').write_text(READ)
        joined = reads_input(self.reads, self.out, False)
        self.assertEqual(joined, self.out / '00_prep/hifi_reads.fastq')
        self.assertEqual(joined.read_text(), READ * 2)

    def test_file_with_other_suffix_is_rejected(self):
        other = self.dir / 'notes.txt'
        other.write_text('x\n')
        with self.assertRaises(UsageError) as caught:
            reads_path_input(other)
        self.assertIn('not a HiFi FASTQ', str(caught.exception))

    def test_empty_file_is_rejected(self):
        empty = self.dir / 'empty.fastq'
        empty.touch()
        with self.assertRaises(UsageError) as caught:
            reads_path_input(empty)
        self.assertIn('empty', str(caught.exception))

    def test_missing_path_is_rejected(self):
        with self.assertRaises(UsageError) as caught:
            reads_path_input(self.dir / 'nope.fastq')
        self.assertIn('not found', str(caught.exception))

    def test_mixed_compression_is_rejected(self):
        (self.reads / 'b.fastq.gz').write_bytes(b'\x1f\x8b\x08\x00\x00\x00\x00\x00')
        with self.assertRaises(UsageError) as caught:
            reads_input(self.reads, self.out, True)
        self.assertIn('mixes gzipped and uncompressed', str(caught.exception))


if __name__ == '__main__':
    unittest.main()
