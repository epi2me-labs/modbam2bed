import itertools
import unittest

NOT_G = 'ACTMWYH'
NOT_C = 'AGTRWKD'
MAYBE_G = 'GRSKVDBN'
MAYBE_C = 'CMQYVHBN'


from modbampy import libbam

class MotifTest(unittest.TestCase):

    def test_001_is_cpg(self):
        is_cpg_fwd = libbam.is_cpg_fwd
        is_cpg_rev = libbam.is_cpg_rev

        # forward
        assert is_cpg_fwd(3, 8, b"AAACGAAA") == True
        assert is_cpg_fwd(3, 8, b"AAACTAAA") == False
        # reverse
        assert is_cpg_rev(4, 8, b"AAACGAAA") == True
        assert is_cpg_rev(4, 8, b"AAACTAAA") == False
        # end, but overrun
        assert is_cpg_fwd(6, 7, b"AAAAAACG") == False
        assert is_cpg_fwd(6, 8, b"AAAAAACG") == True
        # don't break
        assert is_cpg_rev(0, 8, b"AAAAAACG") == False


    def test_010_is_chh(self):
        is_chh_fwd = libbam.is_chh_fwd
        is_chh_rev = libbam.is_chh_rev

        # forward
        for b1, b2 in itertools.product(NOT_G, repeat=2):
            assert is_chh_fwd(3, 8, f"AAAC{b1}{b2}AA".encode()) == True
        for b1, b2 in itertools.product(MAYBE_G, repeat=2):
            assert is_chh_fwd(3, 8, f"AAAC{b1}{b2}AA".encode()) == False
        for b1, b2 in itertools.product(NOT_G, MAYBE_G):
            assert is_chh_fwd(3, 8, f"AAAC{b1}{b2}AA".encode()) == False
            assert is_chh_fwd(3, 8, f"AAAC{b2}{b1}AA".encode()) == False

        # reverse
        for b1, b2 in itertools.product(NOT_C, repeat=2):
            assert is_chh_rev(5, 8, f"AAA{b1}{b2}GAA".encode()) == True
        for b1, b2 in itertools.product(MAYBE_C, repeat=2):
            assert is_chh_rev(5, 8, f"AAA{b1}{b2}GAA".encode()) == False
        for b1, b2 in itertools.product(NOT_C, MAYBE_C):
            assert is_chh_rev(5, 8, f"AAA{b1}{b2}GAA".encode()) == False
            assert is_chh_rev(5, 8, f"AAA{b2}{b1}GAA".encode()) == False

        # end, but overrun
        assert is_chh_fwd(5, 7, b"AAAAACHH") == False
        assert is_chh_fwd(5, 8, b"AAAAACHH") == True

        # don't break
        for i in 0, 1:
            assert is_chh_rev(5, 7, b"AAAAACHH") == False


    def test_020_is_chg(self):
        is_chg_fwd = libbam.is_chg_fwd
        is_chg_rev = libbam.is_chg_rev

        # forward
        for b1 in NOT_G:
            assert is_chg_fwd(3, 8, f"AAAC{b1}GAA".encode()) == True
        for b1 in MAYBE_G:
            assert libbam.is_chg_fwd(3, 8, f"AAAC{b1}GAA".encode()) == False
        for b1, b2 in itertools.product(MAYBE_G, NOT_G + MAYBE_G):
            assert is_chg_fwd(3, 8, f"AAAC{b1}{b2}AA".encode()) == False

        # reverse
        for b1 in NOT_C:
            assert is_chg_rev(5, 8, f"AAAC{b1}GAA".encode()) == True
        for b1 in MAYBE_C:
            assert is_chg_rev(5, 8, f"AAAC{b1}GAA".encode()) == False
        for b1, b2 in itertools.product(NOT_C + MAYBE_C, MAYBE_C):
            assert is_chg_rev(5, 8, f"AAA{b1}{b2}GAA".encode()) == False

        # end, but overrun
        assert is_chg_fwd(5, 7, b"AAAAACHG") == False
        assert is_chg_fwd(5, 8, b"AAAAACHG") == True

        # don't break
        for i in 0, 1:
            assert is_chg_rev(5, 7, b"AAAAACHG") == False
