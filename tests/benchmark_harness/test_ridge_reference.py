import copy
import unittest
from ridge_reference import ridge_fixtures, assert_ridge_result

class RidgeOracleTests(unittest.TestCase):
    def test_analytic_solutions_are_stationary(self):
        for case in ridge_fixtures():
            if not case['invalid']:
                with self.subTest(case=case['name']): assert_ridge_result(case['expected'],case)

    def test_zero_cannot_pass_for_tiny_nonzero_answer(self):
        case=next(c for c in ridge_fixtures() if c['expected'] and c['expected']['values'] and 0<abs(c['expected']['values'][0])<1e-200)
        wrong=copy.deepcopy(case['expected']);wrong['values'][0]=0
        with self.assertRaises(AssertionError): assert_ridge_result(wrong,case)

    def test_unregularized_answer_is_rejected(self):
        case=ridge_fixtures()[0]
        with self.assertRaises(AssertionError): assert_ridge_result(dict(rows=1,cols=1,values=[1.5]),case)
