import json,subprocess,sys,unittest,random
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'benchmarks'))
import sparse
from cholesky_reference import fixtures,analyze,factor,check_result
from sparse_reference import dense_csr,multiply,protocol
class CholeskyTests(unittest.TestCase):
 def run_case(self,case,iterations=0):
  args,stdin=protocol(case,iterations)
  return subprocess.run([sys.executable,str(ROOT/'ports/python/runner.py'),*args],input=stdin,text=True,capture_output=True)
 def test_shared_fixtures(self):
  for case in fixtures():
   with self.subTest(case=case['name']):
    r=self.run_case(case)
    if case['invalid']:self.assertNotEqual(r.returncode,0)
    else:self.assertEqual(r.returncode,0,r.stderr);check_result(json.loads(r.stdout),case)
 def test_random_spd_reconstruction(self):
  rng=random.Random(2026)
  for n in range(1,13):
   edges=[[0.]*n for _ in range(n)]
   for i in range(n):
    for j in range(i):
     if rng.random()<.25:edges[i][j]=edges[j][i]=rng.uniform(-1,1)
   for i in range(n):edges[i][i]=1+sum(abs(v) for v in edges[i])
   a=dense_csr(edges);case=dict(op='chol_factor',a=a,b=[0]*n,expected=sum(factor(a),[]))
   r=self.run_case(case);self.assertEqual(r.returncode,0,r.stderr);check_result(json.loads(r.stdout),case)
 def test_workload_checksums(self):
  from cholesky_bench import cases
  for problem in ['grid','scrambled']:
   for _,_,case,stats in cases(3,problem):
    self.assertGreater(stats['factor_nnz'],0)
    r=self.run_case(case);self.assertEqual(r.returncode,0,r.stderr);check_result(json.loads(r.stdout),case)
    r=self.run_case(case,2);self.assertEqual(r.returncode,0,r.stderr);self.assertAlmostEqual(json.loads(r.stdout)['checksum'],2*sum(case['expected']),delta=1e-7)
 def test_structure_and_cancellation(self):
  star=dense_csr([[5,1,1,1],[1,3,0,0],[1,0,3,0],[1,0,0,3]])
  self.assertEqual(sum(k>=0 for k in analyze(star)[2]),3)
  self.assertEqual([k for k in analyze(star)[2] if k>=0],[0,0,0])
  a=dense_csr([[4,2,2],[2,5,1],[2,1,5]])
  self.assertEqual(factor(a)[2][-2],0)
 def test_public_fields(self):
  from publication import PublicSanitizer
  data=dict(methodology={'storage':'Logical factor bytes only'},implementations=[],results=[dict(factor_nnz=10,fill_count=3,logical_factor_bytes=200,hostname='private')])
  clean=PublicSanitizer().report(data);row=clean['results'][0]
  self.assertEqual(clean['methodology']['storage'],'Logical factor bytes only')
  self.assertEqual(row['factor_nnz'],10);self.assertEqual(row['fill_count'],3);self.assertNotIn('hostname',row)
 def test_refresh_preserves_saved_bundle_and_measurements(self):
  import tempfile
  from cholesky_bench import render
  from refresh_reports import embedded
  data={'implementations':[],'results':[dict(implementation='c',operation='factor',ordering='rcm',problem='grid',size=4,factor_nnz=67,fill_count=27,logical_factor_bytes=1208,median_ns=123.456,status='passed')]}
  live={'available':False,'reason':'Saved offline report'}
  with tempfile.TemporaryDirectory() as directory:
   path=Path(directory)/'index.html';render(data,path,live_override=live)
   self.assertEqual(embedded(path,'live'),live)
   self.assertEqual(embedded(path,'data')['results'],data['results'])
