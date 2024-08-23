import os
import unittest

from modbampy import ModBam

test_bam = os.path.join(os.path.dirname(__file__), "..", "test_data", "400ecoli.bam")
tag_codes_bam = os.path.join(os.path.dirname(__file__), "..", "test_data", "tag_codes.bam")

class MotifTest(unittest.TestCase):

    def test_001_mapping_quality(self):
        t = 0
        i = 0
        with ModBam(test_bam) as bam:
            for r in bam.reads("ecoli1", 0, 4000000):
                t += r.mapping_quality
                i += 1
        assert(t / i == 60.0)

    def test_002_complete_parse(self):
        n_mods = 0
        with ModBam(test_bam) as bam:
            for r in bam.reads("ecoli1", 0, 4000000):
                for site in r.mod_sites:
                    n_mods += 1
        assert n_mods == 3259

    def test_010_tag_parse(self):
        with ModBam(tag_codes_bam) as bam:
            for r in bam.reads("ecoli1", 0, 4000000):
                # read name has expected modified base code
                expected_tag = r.query_name.split("_")[-1]
                try:
                    expected_tag = int(expected_tag)
                except ValueError:
                    pass
                for site in r.mod_sites:
                    print(site)
                    assert site.mbase == expected_tag

    def test_020_pileup(self):
        with ModBam(test_bam) as bam:
            bam.pileup("ecoli1", 105000, 105100)
        with ModBam(tag_codes_bam) as bam:
            with self.assertRaises(ValueError):
                bam.pileup("ecoli1", 0, 4000000, mod_base=27551)
            bam.pileup("ecoli1", 0, 4000000, mod_base=27551, canon_base="C")

    def test_030_pythonic(self):
        with ModBam(test_bam) as bam:
            reads = list(bam.reads("ecoli1", 0, 4000000))
        qnames = set(x.query_name for x in reads)
        assert(len(qnames) > 1)

    def test_040_nonexisting_chrom(self):
        with ModBam(test_bam) as bam:
            reads = list(bam.reads("ecoli1xx", 0, 4000000))
        assert(reads == [])
