"""False-duplication boundary and pair values are Y / N / NA."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'scripts'))
import false_duplication as fd  # noqa: E402


class Pair:
    def __init__(self, duplication_type, orientation='DIRECT'):
        self.duplication_type = duplication_type
        self.orientation = orientation


def statuses(c1l, c1r, c2l, c2r):
    return {'copy1_left': c1l, 'copy1_right': c1r, 'copy2_left': c2l, 'copy2_right': c2r}


class FalseDuplicationCodeTests(unittest.TestCase):
    def test_boundary_status(self):
        ev = fd.BoundaryEvidence()
        self.assertEqual(fd.boundary_status(ev, 2, 2), 'NA')
        ev.left_reaching.update({'a', 'b'})
        ev.right_reaching.update({'a', 'b'})
        self.assertEqual(fd.boundary_status(ev, 2, 2), 'N')
        ev.spanning_reads.add('a')
        self.assertEqual(fd.boundary_status(ev, 2, 2), 'NA')
        ev.spanning_reads.add('b')
        self.assertEqual(fd.boundary_status(ev, 2, 2), 'Y')

    def test_pair_calls(self):
        tandem = Pair('TANDEM')
        self.assertEqual(fd.classify_pair(tandem, statuses('Y', 'N', 'Y', 'Y')),
                         ('Y_false_duplication_pattern', 'Y', 'Y'))
        self.assertEqual(fd.classify_pair(tandem, statuses('NA', 'N', 'Y', 'Y')),
                         ('NA_insufficient_reads', 'N', 'NA'))
        self.assertEqual(fd.classify_pair(tandem, statuses('NA', 'Y', 'NA', 'NA')),
                         ('N_tandem_junction_spanned', 'Y', 'N'))
        other = Pair('INTERSPERSED')
        self.assertEqual(fd.classify_pair(other, statuses('Y', 'Y', 'N', 'N')),
                         ('Y_false_duplication_pattern', 'Y', 'Y'))
        self.assertEqual(fd.classify_pair(other, statuses('Y', 'Y', 'Y', 'Y')),
                         ('N_no_false_duplication_pattern', 'Y', 'N'))

    def test_defining_boundaries(self):
        row = {'LR_FALSE_DUPLICATION': 'Y', 'DUPLICATION_TYPE': 'INTERSPERSED', 'pair_id': 'p1',
               'LR_COPY1_LEFT_SUPPORTED': 'Y', 'LR_COPY1_RIGHT_SUPPORTED': 'Y',
               'LR_COPY2_LEFT_SUPPORTED': 'N', 'LR_COPY2_RIGHT_SUPPORTED': 'N'}
        self.assertEqual(fd.defining_boundaries(row), {'copy2_left', 'copy2_right'})


if __name__ == '__main__':
    unittest.main()
