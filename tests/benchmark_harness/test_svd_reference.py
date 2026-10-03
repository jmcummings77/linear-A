import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"benchmarks"))
import copy
import unittest
from svd_reference import assert_svd_result
from reference import matrix

class SvdOracleTests(unittest.TestCase):
    def test_sign_changes_are_valid_but_bad_factors_and_spectra_are_not(self):
        a=matrix(2,2,[3,0,0,2])
        good=matrix(5,2,[1,0,0,1,3,2,1,0,0,1])
        assert_svd_result(good,a,[3,2])
        changed=copy.deepcopy(good);changed['values'][0]=-1;changed['values'][6]=-1
        assert_svd_result(changed,a,[3,2])
        for index,value in [(0,2),(4,-3),(4,4),(5,4),(6,-1),(0,float('nan'))]:
            bad=copy.deepcopy(good);bad['values'][index]=value
            with self.assertRaises(AssertionError):assert_svd_result(bad,a,[3,2])
