import os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vre.scoring import score_row, tier_for
from vre.cli import run

D = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sample_data")
FILES = [os.path.join(D, f) for f in sorted(os.listdir(D))]


class T(unittest.TestCase):
    def test_pathogenic_lof_rare_is_very_high(self):
        r = score_row({"gene": "G", "clinvar": "Pathogenic", "consequence": "stop_gained",
                       "gnomad_af": "0.00001", "cadd": "38", "matched_hpo": "HP:1,HP:2"})
        self.assertEqual(r.tier, "Very High")

    def test_benign_common_is_very_low(self):
        r = score_row({"clinvar": "Benign", "consequence": "synonymous_variant", "gnomad_af": "0.2"})
        self.assertEqual(r.tier, "Very Low")

    def test_likely_pathogenic_not_scored_as_pathogenic(self):
        self.assertTrue(any("+32" in x for x in score_row({"clinvar": "Likely_pathogenic"}).reasons))

    def test_missing_columns_do_not_crash(self):
        self.assertIn(score_row({}).tier, {"Very High", "High", "Medium", "Low", "Very Low"})

    def test_every_variant_has_reasons(self):
        for r in run(FILES):
            self.assertGreaterEqual(len(r.reasons), 4)

    def test_highest_first(self):
        res = run(FILES)
        self.assertEqual(res[0].tier, "Very High")

    def test_tier_cutoffs(self):
        self.assertEqual((tier_for(75), tier_for(74), tier_for(-5)), ("Very High", "High", "Very Low"))


if __name__ == "__main__":
    unittest.main()
