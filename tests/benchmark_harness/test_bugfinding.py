"""Test generators and reduction with deliberate faulty runners, not only good outputs."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'benchmarks'))
import bugfinding as hunt

class BugFindingTests(unittest.TestCase):
    def test_seed_reproduces_families_and_independent_oracles(self):
        a=list(hunt.generated(73));self.assertEqual(a,list(hunt.generated(73)))
        self.assertNotEqual(a,list(hunt.generated(74)))
        self.assertEqual(len(a),10)
        for seed in range(12):
            for case in hunt.generated(seed):hunt.oracle(case)

    def test_bad_protocol_nonfinite_and_timeout_cannot_pass(self):
        case=next(hunt.generated(4));implementation=dict(runner=['runner'],env={})
        for output in ('not JSON','null','{"rows":true,"cols":1,"values":[0]}'):
            with patch.object(hunt.harness,'execute',return_value=SimpleNamespace(returncode=0,stdout=output,stderr='')):
                with self.assertRaises(hunt.Failure):hunt.evaluate(implementation,case)
        with patch.object(hunt.harness,'execute',side_effect=subprocess.TimeoutExpired('runner',30)):
            with self.assertRaises(hunt.Failure) as error:hunt.evaluate(implementation,case)
            self.assertEqual(error.exception.kind,'timeout')

    def test_reducer_finds_small_counterexample_and_preserves_original(self):
        case=dict(generator_version=1,seed=7,family='injected multiply defect',op='multiply',
                  a=hunt.matrix(3,3,[4]*9),b=hunt.matrix(3,3,[2]*9))
        original=copy.deepcopy(case)
        def fails(c):return any(v!=0 for v in hunt.oracle(c)['values'])
        reduced,attempts=hunt.reduce_case(case,fails,200)
        self.assertEqual(case,original);self.assertTrue(fails(reduced));self.assertLessEqual(attempts,200)
        self.assertEqual((reduced['a']['rows'],reduced['a']['cols'],reduced['b']['cols']),(1,1,1))
        self.assertEqual(reduced['a']['values'],[1]);self.assertEqual(reduced['b']['values'],[1])
        unchanged,attempts=hunt.reduce_case(case,fails,0)
        self.assertEqual(unchanged,case);self.assertEqual(attempts,0)

    def test_invalid_singular_and_nonsymmetric_candidates_never_reach_predicate(self):
        case=next(c for c in hunt.generated(73) if c['op']=='solve_cholesky')
        seen=[]
        def predicate(candidate):hunt.oracle(candidate);seen.append(candidate);return False
        hunt.reduce_case(case,predicate,50)
        self.assertTrue(seen)
        invalid=copy.deepcopy(case);invalid['a']['values']=[0]*len(invalid['a']['values'])
        with self.assertRaises(ValueError):hunt.oracle(invalid)

    def test_saved_failure_is_sanitized_replayable_and_does_not_claim_global_minimum(self):
        case=dict(generator_version=1,seed=9,family='injected transpose defect',op='transpose_twice',
                  a=hunt.matrix(2,2,[1,2,3,4]),b=None)
        def broken(implementation,c):
            if c['a']['values']:raise hunt.Failure('numerical','bad result at /Users/private-name/source.c')
            return hunt.oracle(c)
        with tempfile.TemporaryDirectory() as directory,patch.object(hunt,'evaluate',side_effect=broken):
            checks,failures=hunt.run_cases([dict(id='c')],[case],Path(directory),25)
            artifact=json.loads(next(Path(directory).glob('failure-*.json')).read_text())
        self.assertFalse(checks[0]['passed']);self.assertEqual(len(failures),1)
        self.assertNotIn('/Users/',json.dumps(artifact));self.assertFalse(artifact['globally_minimal'])
        self.assertTrue(artifact['reproduced_after_reduction']);hunt.oracle(artifact['reduced'])

    def test_replay_rejects_unbounded_or_malformed_inputs(self):
        case=next(hunt.generated(4))
        for update in [dict(rows=1000000),dict(rows=True),dict(values=[float('nan')]),dict(values='oops')]:
            bad=copy.deepcopy(case);bad['a'].update(update)
            with self.assertRaises(ValueError):hunt.validate(bad)

    def test_property_checks_reject_compensating_transpose_defects(self):
        case=dict(generator_version=1,seed=1,family='transpose property',op='transpose_twice',
                  a=hunt.matrix(2,2,[1,2,3,4]),b=None)
        # Returning the input twice satisfies involution but is not transpose.
        with patch.object(hunt,'request',return_value=case['a']):
            with self.assertRaises(hunt.Failure) as failure:hunt.evaluate({},case)
        self.assertEqual(failure.exception.kind,'numerical')

    def test_differential_failure_rechecks_both_ports_during_reduction(self):
        case=dict(generator_version=1,seed=1,family='differential path',op='determinant',
                  a=hunt.matrix(2,2,[1,0,0,1]),b=None)
        calls=[]
        def differing(implementation,c):
            calls.append(implementation['id'])
            return {'value':1 if implementation['id']=='c' else 2}
        with tempfile.TemporaryDirectory() as directory,patch.object(hunt,'evaluate',side_effect=differing):
            checks,failures=hunt.run_cases([{'id':'c'},{'id':'python'}],[case],Path(directory),10)
        self.assertEqual(failures[0]['kind'],'differential')
        self.assertTrue(failures[0]['reproduced_after_reduction'])
        self.assertGreater(calls.count('c'),1);self.assertGreater(calls.count('python'),1)
