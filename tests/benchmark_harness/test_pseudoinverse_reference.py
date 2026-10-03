import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'benchmarks'))
from pseudoinverse_reference import pseudoinverse_fixtures, assert_inverse_result

class PseudoinverseOracleTests(unittest.TestCase):
    def test_analytic_answers_and_invariants_agree(self):
        for case in pseudoinverse_fixtures():
            if not case['invalid']:
                with self.subTest(case=case['name']): assert_inverse_result(case['expected'],case)

    def test_bad_answers_are_rejected_including_tiny_nonzeros(self):
        for case in pseudoinverse_fixtures():
            if case['invalid']: continue
            result=copy.deepcopy(case['expected'])
            index=next((i for i,x in enumerate(result['values']) if x),None)
            if index is not None:
                result['values'][index]=0
                with self.subTest(case=case['name']),self.assertRaises(AssertionError): assert_inverse_result(result,case)
