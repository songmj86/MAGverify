"""High-confidence clipping requires validated clipping and the stringent thresholds.

Criterion (Supplementary Methods): a candidate site from
anvi-script-find-misassemblies (clipping ratio > 0.8) is high-confidence only when
it is HiFi-supported (>= 2 distinct clipped reads, no spanning read) AND has
anvi'o coverage >= 10 with a clipping ratio of 1.0. The threshold pair alone
(Clipping_candidate=Y) is not sufficient.
"""
from pathlib import Path
import csv
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))
import clipping_summarize as cs  # noqa: E402

FIELDS = sorted({
    "contig", "length", "pos_1based", "coverage", "clipping_count",
    "clipping_ratio", "default_candidate", "Clipping_candidate", "Clipping_ge1",
    "Clipping_analyzable", "Clipping", "High_confidence_clipping", "region_type",
    "n_unique_all_clipped_qnames", "n_unique_exact_spanning_qnames",
    "exact_spanning_min_flank_bp_min", "exact_spanning_min_flank_bp_median",
    "exact_spanning_min_flank_bp_max",
    "n_unique_spanning_with_I_qnames", "spanning_I_length_bp_min",
    "spanning_I_length_bp_median", "spanning_I_length_bp_max",
    "n_unique_spanning_with_D_qnames", "spanning_D_length_bp_min",
    "spanning_D_length_bp_median", "spanning_D_length_bp_max",
    "n_unique_spanning_with_N_qnames", "spanning_N_length_bp_min",
    "spanning_N_length_bp_median", "spanning_N_length_bp_max",
    "n_unique_spanning_with_any_IDN_qnames",
})


def site(**overrides):
    row = {field: 'NA' for field in FIELDS}
    row.update(
        contig='contig_1', length=100000, pos_1based=500,
        coverage=20, clipping_count=20, clipping_ratio='1.000000',
        default_candidate='Y', Clipping_candidate='Y',
        Clipping_ge1='Y', Clipping_analyzable='Y',
        Clipping='Y', High_confidence_clipping='Y',
        region_type='INTERNAL', n_unique_all_clipped_qnames=20,
    )
    for field in FIELDS:
        if field.startswith('n_unique_') and field != 'n_unique_all_clipped_qnames':
            row[field] = 0
    row.update(overrides)
    return row


class HighConfidenceClippingCriterion(unittest.TestCase):
    def load(self, rows):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'clipping.site.tsv'
            with open(path, 'w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t')
                writer.writeheader()
                writer.writerows(rows)
            return cs.load_clipping(str(path), {'contig_1': 100000})

    def test_validated_site_meeting_thresholds_is_accepted(self):
        _fields, records = self.load([site()])
        self.assertEqual(records[0]['High_confidence_clipping'], 'Y')
        self.assertEqual(records[0]['Clipping'], 'Y')

    def test_thresholds_without_validated_clipping_are_rejected(self):
        # coverage 20x, ratio 1.0 (Clipping_candidate=Y) but not HiFi-supported
        with self.assertRaises(SystemExit) as caught:
            self.load([site(Clipping='N', Clipping_ge1='N', High_confidence_clipping='Y')])
        self.assertIn('High_confidence_clipping=Y requires Clipping=Y', str(caught.exception))

    def test_candidate_may_meet_thresholds_without_being_high_confidence(self):
        # the threshold pair alone leaves High_confidence_clipping=N
        _fields, records = self.load([
            site(Clipping='N', Clipping_ge1='N', Clipping_candidate='Y', High_confidence_clipping='N')
        ])
        self.assertEqual(records[0]['Clipping_candidate'], 'Y')
        self.assertEqual(records[0]['High_confidence_clipping'], 'N')

    def test_unanalyzable_site_is_na(self):
        with self.assertRaises(SystemExit) as caught:
            self.load([site(Clipping_analyzable='N', Clipping='NA', High_confidence_clipping='N')])
        self.assertIn('requires High_confidence_clipping=NA', str(caught.exception))


if __name__ == '__main__':
    unittest.main()
