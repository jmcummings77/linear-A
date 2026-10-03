import copy,json,subprocess,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'benchmarks'))
from ilu_reference import fixtures,reference_apply
from sparse_reference import protocol,check_result as check_simple
from gmres_reference import check_result as check_gmres,dense_solve
from publication import PublicSanitizer
class IluTests(unittest.TestCase):
    def test_python_against_independent_exact_oracle(self):
        for case in fixtures():
            with self.subTest(case=case['name']):
                args,stdin=protocol(case);r=subprocess.run([sys.executable,str(ROOT/'ports/python/runner.py'),*args],input=stdin,text=True,capture_output=True)
                if case['invalid']:self.assertNotEqual(r.returncode,0)
                else:
                    self.assertEqual(r.returncode,0,r.stderr)
                    (check_gmres if case['op']=='gmres' else check_simple)(json.loads(r.stdout),case)
    def test_dropped_fill_is_not_silently_dense_lu(self):
        case=fixtures()[0];self.assertNotEqual(reference_apply(case['a'],case['b']),dense_solve(case['a'],case['b']))
        explicit=fixtures()[3];actual=reference_apply(explicit['a'],explicit['b']);expected=dense_solve(explicit['a'],explicit['b'])
        for x,y in zip(actual,expected):self.assertAlmostEqual(x,y,places=14)
    def test_factor_storage_is_public_but_environment_is_not(self):
        row=PublicSanitizer().report({'results':[{'logical_preconditioner_bytes':1234,'environment':{'TOKEN':'secret'}}]})['results'][0]
        self.assertEqual(row['logical_preconditioner_bytes'],1234);self.assertNotIn('environment',row)
