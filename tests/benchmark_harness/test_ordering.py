import json,subprocess,sys,unittest,random
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'benchmarks'))
from ordering_reference import fixtures,permute,rcm,bandwidth
from sparse_reference import protocol,dense_csr,multiply
import sparse
class OrderingTests(unittest.TestCase):
 def test_shared_fixtures(self):
  for case in fixtures():
   with self.subTest(case=case['name']):
    args,stdin=protocol(case);r=subprocess.run([sys.executable,str(ROOT/'ports/python/runner.py'),*args],input=stdin,text=True,capture_output=True)
    if case['invalid']:self.assertNotEqual(r.returncode,0)
    else:self.assertEqual(r.returncode,0,r.stderr);sparse.check_result(json.loads(r.stdout),case)
 def test_graph_convention_and_determinism(self):
  self.assertEqual(rcm(dense_csr([[1,1,0],[0,1,1],[0,0,1]])),[2,1,0])
  self.assertEqual(rcm(dense_csr([[1,0,0],[0,2,0],[0,0,3]])),[2,1,0])
 def test_random_permutation_similarity_and_inverse(self):
  import importlib.util
  spec=importlib.util.spec_from_file_location('ordering_python_sparse',ROOT/'ports/python/sparse.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);CSRMatrix=module.CSRMatrix
  rng=random.Random(2026)
  for n in range(1,14):
   a=dense_csr([[rng.randint(-8,8) if rng.random()<.2 or i==j else 0 for j in range(n)] for i in range(n)])
   p=list(range(n));rng.shuffle(p);q=permute(a,p)
   obj=CSRMatrix(n,n,a['offsets'],a['indices'],a['values']);actual=obj.permute_symmetric(p)
   self.assertEqual(actual.values,q['values']);self.assertEqual(actual.column_indices,q['indices']);self.assertEqual(obj.reverse_cuthill_mckee(),rcm(a))
   x=[rng.random() for _ in range(n)];y=CSRMatrix.permute_vector(p,x)
   self.assertEqual(CSRMatrix.permute_vector(p,y,True),x)
   z=CSRMatrix.permute_vector(p,actual.matvec(y),True)
   for v,w in zip(z,obj.matvec(x)):self.assertAlmostEqual(v,w,places=12)
   inv=[p.index(i) for i in range(n)];back=actual.permute_symmetric(inv)
   self.assertEqual(back.values,obj.values);self.assertEqual(back.column_indices,obj.column_indices)
 def test_every_benchmark_workload_uses_valid_options_and_original_residual(self):
  from ordering_bench import cases
  for problem in ('grid','scrambled'):
   for _,_,case,_ in cases(4,problem):
    with self.subTest(problem=problem,op=case['op']):
     args,stdin=protocol(case);r=subprocess.run([sys.executable,str(ROOT/'ports/python/runner.py'),*args],input=stdin,text=True,capture_output=True)
     self.assertEqual(r.returncode,0,r.stderr);sparse.check_result(json.loads(r.stdout),case)
 def test_benchmark_calibration_and_samples_include_iteration_checksum(self):
  import contextlib,io
  import run,ordering_bench
  with contextlib.redirect_stdout(io.StringIO()):rows=ordering_bench.benchmark([run.build_one('python',no_build=True)],[3],1)
  self.assertEqual(len(rows),16)
  self.assertTrue(all(r['status']=='passed' for r in rows),[r.get('error') for r in rows if r['status']!='passed'])
  self.assertTrue(all(len(r['samples'])==1 for r in rows))
 def test_public_ordering_metadata_has_no_environment_or_paths(self):
  from publication import PublicSanitizer
  data=PublicSanitizer().report({'results':[{'implementation':'c','ordering':'rcm','problem':'grid','bandwidth':12,'environment':{'TOKEN':'secret'}}]})
  self.assertEqual(data['results'][0],{'implementation':'c','ordering':'rcm','problem':'grid','bandwidth':12})
 def test_report_refresh_preserves_embedded_live_bundle_and_measurements(self):
  import tempfile
  from ordering_bench import render
  from refresh_reports import embedded
  data={'implementations':[],'results':[{'implementation':'c','operation':'total','ordering':'rcm','problem':'grid','size':12,'bandwidth':12,'median_ns':123,'status':'passed'}]}
  live={'available':False,'reason':'No verified build'}
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'index.html';render(data,path,live_override=live)
   self.assertEqual(embedded(path,'live'),live);self.assertEqual(embedded(path,'data')['results'],data['results'])
