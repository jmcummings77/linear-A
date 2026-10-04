import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'benchmarks'))
from multigrid_reference import poisson,transfer,reference,fixtures
from multigrid_bench import cases,stencil_ic0
from ic0_reference import factor
from publication import PublicSanitizer

class MultigridReferenceTests(unittest.TestCase):
 def test_transfer_boundary_and_adjoint(self):
  p=transfer(3)
  self.assertEqual([float(row[0]) for row in p],[.25,.5,.25,.5,1,.5,.25,.5,.25])
  # Bilinear prolongation and its transpose, with homogeneous outside boundaries.
  p=transfer(7);x=list(range(9));y=list(range(49))
  px=[sum(a*b for a,b in zip(row,x)) for row in p];pty=[sum(p[i][j]*y[i] for i in range(49)) for j in range(9)]
  self.assertEqual(sum(a*b for a,b in zip(px,y)),sum(a*b for a,b in zip(x,pty)))
 def test_stencil_factor_reference(self):
  for w in (1,3,7):
   actual=stencil_ic0(w);expected=factor(poisson(w),exact=False)
   self.assertEqual(actual[:2],expected[:2])
   for a,b in zip(actual[2],expected[2]):self.assertAlmostEqual(a,b,places=13)
 def test_fixture_expected_values_are_finite_and_setup_is_distinct(self):
  import math
  self.assertEqual(len(fixtures()),21)
  for f in fixtures():
   if f['expected'] is not None:self.assertTrue(all(math.isfinite(v) for v in f['expected']))
  rows=list(cases(7));byop={r[0]:r for r in rows}
  self.assertEqual(byop['mg_setup'][2]['expected'],[49,3])
  self.assertEqual(byop['mg_solve'][2]['options']['jacobi'],4)
  self.assertEqual(byop['mg_total'][2]['options']['jacobi'],5)
  self.assertEqual(byop['mg_solve'][3]['logical_factor_bytes'],8)
  self.assertGreater(byop['mg_solve'][3]['logical_workspace_bytes'],49*8)

class MultigridReportTests(unittest.TestCase):
 def test_refresh_preserves_measurements_and_embedded_runtime(self):
  import json,tempfile
  from unittest.mock import patch
  from multigrid_bench import render
  from refresh_reports import embedded
  root=Path(__file__).resolve().parents[2]
  source=root/'benchmarks/reports/multigrid/index.html'
  data=embedded(source,'data');live=embedded(source,'live')
  before=json.dumps(data,sort_keys=True);runtime=json.dumps(live,sort_keys=True)
  with tempfile.TemporaryDirectory() as d,patch('sparse.live_bundle',side_effect=AssertionError('refresh must not build a runtime')):
   target=Path(d)/'index.html';render(data,target,live_override=live,published=True)
   self.assertEqual(embedded(target,'data'),data)
   self.assertEqual(embedded(target,'live'),live)
  self.assertEqual(json.dumps(data,sort_keys=True),before)
  self.assertEqual(json.dumps(live,sort_keys=True),runtime)
  self.assertIn('multigrid',live['capabilities'])
  self.assertEqual(len(live['fixtures']),277)

class MultigridCaptureTests(unittest.TestCase):
 def test_clean_capture_supports_multigrid_without_a_nonexistent_flag(self):
  from reproduce import checked_arguments
  from check_provenance import RUNNERS
  self.assertIn("benchmarks/multigrid_bench.py",RUNNERS)
  self.assertEqual(checked_arguments("multigrid_bench.py",["--sizes","7"]),["--sizes","7"])
