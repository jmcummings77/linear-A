import unittest,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"benchmarks"))
from ic0_reference import fixtures
class IC0HarnessTests(unittest.TestCase):
 def test_shared_fixtures(self):
  import subprocess,json
  from sparse_reference import protocol
  from sparse_reference import check_result
  for case in fixtures():
   with self.subTest(case=case['name']):
    args,data=protocol(case);r=subprocess.run([sys.executable,str(ROOT/'ports/python/runner.py'),*args],input=data,text=True,capture_output=True)
    if case['invalid']:self.assertNotEqual(r.returncode,0)
    else:self.assertEqual(r.returncode,0,r.stderr);check_result(json.loads(r.stdout),case)
 def test_workload_checksums_and_orderings(self):
  import subprocess,json
  from ic0_bench import cases
  from sparse_reference import protocol
  import sparse
  for _,_,case,_ in cases(3):
   args,data=protocol(case);r=subprocess.run([sys.executable,str(ROOT/'ports/python/runner.py'),*args],input=data,text=True,capture_output=True)
   self.assertEqual(r.returncode,0,r.stderr);sparse.check_result(json.loads(r.stdout),case)
   args,data=protocol(case,2);r=subprocess.run([sys.executable,str(ROOT/'ports/python/runner.py'),*args],input=data,text=True,capture_output=True)
   self.assertEqual(r.returncode,0,r.stderr);self.assertAlmostEqual(json.loads(r.stdout)['checksum'],2*sum(case['expected']),delta=1e-7)
 def test_exact_oracle_positive_definite_breakdown(self):
  from ic0_reference import reference,BREAKDOWN
  from sparse_reference import dense_csr
  with self.assertRaises(ValueError):reference(dense_csr(BREAKDOWN))
