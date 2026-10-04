import importlib.util,random,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'benchmarks'))
from amd_reference import amd,fixtures,matrix
from sparse_reference import dense_csr
spec=importlib.util.spec_from_file_location('amd_python_sparse',ROOT/'ports/python/sparse.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
class AMDTests(unittest.TestCase):
 def test_random_bounds_and_deterministic_ports(self):
  rng=random.Random(2026);strict=False
  for n in range(1,25):
   for _ in range(8):
    a=dense_csr([[1 if i==j or rng.random()<.1 else 0 for j in range(n)] for i in range(n)])
    expected,frames=amd(a,True);obj=module.CSRMatrix(n,n,a['offsets'],a['indices'],a['values'])
    self.assertEqual(obj.approximate_minimum_degree(),expected)
    self.assertEqual(obj.approximate_minimum_degree(),expected)
    self.assertEqual(sorted(expected),list(range(n)))
    strict |= any(any(f['estimates'][i]>f['exact'][i] for i in f['exact']) for f in frames)
  self.assertTrue(strict,'exercise approximate rather than exact degree')
 def test_fixture_distinguishes_approximate_from_exact_degree(self):
  c=next(c for c in fixtures() if c['name']=='strict bound changes pivot')
  p,frames=amd(c['a'],True)
  self.assertTrue(any(f['pivot']!=min(f['exact'],key=lambda i:(f['exact'][i],i)) for f in frames))
  self.assertEqual(p,[13,5,0,1,3,7,9,10,16,4,6,2,8,11,12,14,15,17,18])
 def test_families_and_fill(self):
  from cholesky_reference import analyze
  from ordering_reference import permute
  for family in ('grid','scrambled','tree','irregular'):
   a=matrix(5,family);p,frames=amd(a,True);rp,ci,steps=analyze(permute(a,p))
   self.assertEqual(sum(len(f['fill']) for f in frames),sum(k>=0 for k in steps))
  a=matrix(5,'tree');self.assertEqual(sum(len(f['fill']) for f in amd(a,True)[1]),0)
 def test_fixtures(self):
  import sparse,run
  impl=run.build_one('python',no_build=True);sparse.verify(impl,fixtures())
  self.assertEqual(impl['status'],'passed',impl.get('checks'))
 def test_benchmark_cases_and_published_reference_fields(self):
  import json,subprocess,sparse
  from amd_bench import cases,render
  from sparse_reference import protocol
  from publication import PublicSanitizer
  from refresh_reports import embedded
  import tempfile
  for family in ('grid','scrambled','tree','irregular'):
   for operation,ordering,case,stats in cases(2,family):
    args,stdin=protocol(case,2)
    r=subprocess.run([sys.executable,str(ROOT/'ports/python/runner.py'),*args],input=stdin,text=True,capture_output=True)
    self.assertEqual(r.returncode,0,r.stderr)
    self.assertAlmostEqual(json.loads(r.stdout)['checksum'],2*sum(case['expected']),delta=1e-7)
  row=dict(implementation='suitesparse',operation='ordering',ordering='suitesparse',problem='tree',size=8,fill_count=0,factor_nnz=127,logical_factor_bytes=2552,median_ns=123.456,status='passed',samples=[dict(ns_per_op=123.456,iterations=3,checksum=6048)],environment={'secret':'private'})
  data=dict(implementations=[dict(id='suitesparse',toolchain='AMD 3.3.3',library_path='/home/private/libamd.so')],results=[row])
  clean=PublicSanitizer().report(data)
  self.assertNotIn('environment',clean['results'][0]);self.assertNotIn('library_path',clean['implementations'][0])
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'index.html';live={'available':False,'reason':'offline'};render(data,path,live_override=live)
   self.assertEqual(embedded(path,'live'),live);self.assertEqual(embedded(path,'data')['results'][0]['samples'],row['samples'])
if __name__=='__main__':unittest.main()
