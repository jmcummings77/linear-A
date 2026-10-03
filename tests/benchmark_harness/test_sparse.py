import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'benchmarks'))
from sparse_reference import fixtures,check_result,diffusion,multiply
from publication import PublicSanitizer
spec=importlib.util.spec_from_file_location('sparse_harness',ROOT/'benchmarks/sparse.py')
harness=importlib.util.module_from_spec(spec);spec.loader.exec_module(harness)

class SparseHarnessTests(unittest.TestCase):
    def test_grid_symmetry_and_positive_energy(self):
        a=diffusion(8,2)
        self.assertEqual(len(a['values']),5*64-4*8)
        entries={(i,a['indices'][p]):a['values'][p] for i in range(64) for p in range(a['offsets'][i],a['offsets'][i+1])}
        for (i,j),v in entries.items():self.assertEqual(v,entries[j,i])
        # Strict boundary rows plus connected positive edge weights imply SPD.
        self.assertTrue(all(entries[i,i]>0 for i in range(64)))
        for x in [[1.]*64,[(-1.)**i for i in range(64)]]:self.assertGreater(sum(v*y for v,y in zip(x,multiply(a,x))),0)
    def test_oracle_rejects_false_convergence_and_corrupted_frames(self):
        case=fixtures()[0]
        wrong={'rows':1,'cols':6,'values':[0,0,1,0,0,0]}
        with self.assertRaises(AssertionError):check_result(wrong,case)
        case=copy.deepcopy(case);case['options']['capture']=0;case['expected']=None
        with self.assertRaises(AssertionError):check_result(wrong,case)
    def test_sparse_metadata_is_allowlisted(self):
        data={'results':[dict(implementation='c',operation='cg',unknowns=64,nnz=288,logical_dense_bytes=32768,logical_csr_bytes=5128,cg_iterations=12,env={'SECRET':'x'},hostname='private-host')], 'implementations':[dict(id='c',build_commands=[['/Users/alice/bin/clang','/work/project/ports/c/matrix.c']],env={'SECRET':'x'})]}
        result=PublicSanitizer(root='/work/project',home='/Users/alice').report(data)
        self.assertEqual(result['results'][0]['nnz'],288)
        self.assertNotIn('env',result['results'][0]);self.assertNotIn('hostname',result['results'][0])
        self.assertEqual(result['implementations'][0]['build_commands'][0],['clang','ports/c/matrix.c'])
    def test_render_escapes_metadata_without_rewriting_measurements(self):
        data={'results':[],'implementations':[],'revision':'</script><script>bad()</script>','machine':{}}
        original=copy.deepcopy(data)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'index.html';harness.render(data,p);html=p.read_text()
            self.assertNotIn('</script><script>bad()',html)
            self.assertIn('Copyright © 2026 J.M. Cummings.',html)
        self.assertEqual(data,original)
