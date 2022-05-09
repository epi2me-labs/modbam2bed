import os
import unittest

from modbampy import ModBam

test_bam = os.path.join(os.path.dirname(__file__), "..", "test_data", "400ecoli.bam")

class MotifTest(unittest.TestCase):

    def test_001_mapping_quality(self):
        t = 0
        i = 0
        with ModBam(test_bam) as bam:
            for r in bam.reads("ecoli1", 0, 4000000):
                t += r.mapping_quality
                i += 1
        assert(t / i == 60.0)
