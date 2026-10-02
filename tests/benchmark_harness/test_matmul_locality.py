"""Regression checks for the isolated native multiplication experiment."""
import argparse
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('matmul_locality',ROOT/'experiments/matmul-locality/run.py')
study=importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)

class LocalityTests(unittest.TestCase):
    def test_rectangular_exact_reference(self):
        self.assertEqual(study.reference([1,2,3,4,5,6],[7,8,9,10,11,12],2,3,2),[58,64,139,154])
        self.assertEqual(study.reference([],[],3,0,2),[0]*6)
        self.assertEqual(study.reference([],[],0,3,2),[])

    def result(self,iterations=3,elapsed=100,mode='allocated'):
        return dict(iterations=iterations,elapsed_ns=elapsed,checksum=study.expected_checksum((3,5,7),19)*iterations,
                    full_result_verified=True,output_allocations=iterations if mode=='allocated' else 0,
                    requested_output_bytes=iterations*3*7*8 if mode=='allocated' else 0)

    def test_allocation_and_correctness_claims_are_validated(self):
        for update in [{'checksum':99},{'full_result_verified':False},{'iterations':9},{'elapsed_ns':-1},
                       {'output_allocations':0},{'requested_output_bytes':0}]:
            with self.subTest(update=update),patch.object(study,'execute',return_value=json.dumps(dict(self.result(),**update))):
                with self.assertRaises(AssertionError):study.invoke('runner','ijk',(3,5,7),3,19,'allocated')

    def test_coarse_clock_retries_are_bounded_and_account_for_actual_iterations(self):
        with patch.object(study,'execute',side_effect=[json.dumps(self.result(elapsed=0)),json.dumps(self.result(iterations=48))]) as command:
            data=study.invoke('runner','ijk',(3,5,7),3,19,'allocated')
            self.assertEqual(data['iterations'],48)
            self.assertEqual(command.call_args[0][0][6],48)
        with patch.object(study,'execute',return_value=json.dumps(self.result(elapsed=0))):
            with self.assertRaisesRegex(AssertionError,'bounded'):study.invoke('runner','ijk',(3,5,7),3,19,'allocated')

    def test_public_renderer_preserves_numbers_and_escapes_embedded_markup(self):
        data={'median_ns':12.75,'note':'</script><script>alert(1)</script>','command':study.SANITIZER.command(['/private/tmp/tools/clang',ROOT/'experiments/matmul-locality/src/kernels.c'])}
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory);study.render(data,output);text=(output/'index.html').read_text()
        encoded=text.split('<script id="data" type="application/json">')[1].split('</script>')[0]
        decoded=json.loads(encoded)
        self.assertEqual(decoded['median_ns'],12.75)
        self.assertEqual(decoded['command'][0],'clang')
        self.assertEqual(decoded['command'][1],'experiments/matmul-locality/src/kernels.c')
        self.assertNotIn('</script>',encoded)
        self.assertNotIn('/private/tmp',text)

    def test_only_bounded_positive_timing_shapes(self):
        self.assertEqual(study.shape('96x257x65'),(96,257,65))
        for value in ['1x2','0x2x3','1025x1x1','-1x2x3','AxBxC']:
            with self.assertRaises(argparse.ArgumentTypeError):study.shape(value)
