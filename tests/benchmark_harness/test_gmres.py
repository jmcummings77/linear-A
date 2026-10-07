import copy
import json
import subprocess
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'benchmarks'))
from gmres_reference import fixtures,check_result,dense_solve,transport,unpack
from sparse_reference import multiply,protocol,norm
from publication import PublicSanitizer

class GmresTests(unittest.TestCase):
    def packed(self,case):
        args,stdin=protocol(case)
        result=subprocess.run([sys.executable,str(ROOT/'ports/python/runner.py'),*args],input=stdin,text=True,capture_output=True)
        if result.returncode:raise ValueError(result.stderr)
        return json.loads(result.stdout)
    def test_all_shared_fixtures(self):
        for case in fixtures():
            with self.subTest(case=case['name']):
                if case['invalid']:
                    with self.assertRaises((ValueError,ArithmeticError)):self.packed(case)
                else:check_result(self.packed(case),case)
    def test_oracle_rejects_false_estimates_and_convergence(self):
        case=fixtures()[0];actual=self.packed(case);n=case['a']['rows'];length=actual['values'][2]
        bad=copy.deepcopy(actual);bad['values'][3+n+length]=0
        with self.assertRaises(AssertionError):check_result(bad,case)
        bad=copy.deepcopy(actual);bad['values'][3]=0
        with self.assertRaises(AssertionError):check_result(bad,case)
    def test_shifted_jordan_matches_worked_example_and_right_preconditioning(self):
        cases={case['name']:case for case in fixtures()}
        for jacobi in [0,1]:
            case=cases['GMRES shifted Jordan s=0.25 restart 1 Jacobi %d'%jacobi]
            result=check_result(self.packed(case),case)
            diagonal=8 if jacobi else 1
            self.assertAlmostEqual(result['iterates'][1][0],0,places=12)
            self.assertAlmostEqual(result['iterates'][1][1]*diagonal,4/5,places=12)
            self.assertAlmostEqual(result['x'][0],-32/65,places=12)
            self.assertAlmostEqual(result['x'][1]*diagonal,68/65,places=12)
            self.assertAlmostEqual(result['residuals'][1],1/(5**.5),places=12)
            self.assertAlmostEqual(result['residuals'][2],13**.5/65,places=12)
            self.assertEqual(result['restarts'],[1])
        normal=cases['GMRES same spectrum normal identity']
        jordan=cases['GMRES shifted Jordan s=0.25 restart 2 Jacobi 0']
        self.assertEqual(check_result(self.packed(normal),normal)['iterations'],1)
        self.assertEqual(check_result(self.packed(jordan),jordan)['iterations'],2)
    def test_analytic_oracle_rejects_consistent_but_wrong_intermediate_iterate(self):
        for jacobi in [0,1]:
            case=next(case for case in fixtures() if case['name']=='GMRES shifted Jordan s=1 restart 2 Jacobi %d'%jacobi)
            actual=self.packed(case);n=case['a']['rows'];length=int(actual['values'][2])
            result=unpack(actual,n,True)
            frame=list(result['iterates'][1]);frame[1]*=2
            value=norm([b-y for b,y in zip(case['b'],multiply(case['a'],frame))])
            bad=copy.deepcopy(actual)
            # Keep the actual/estimated histories consistent with the corrupted
            # captured frame. Residual consistency alone cannot reject this.
            bad['values'][3+n+1]=value
            bad['values'][3+n+length+1]=value
            first_frame=len(bad['values'])-n*length+n
            bad['values'][first_frame:first_frame+n]=frame
            with self.assertRaisesRegex(AssertionError,'incorrect analytic iterate'):
                check_result(bad,case)
    def test_analytic_oracle_rejects_wrong_noninitial_residual_estimate(self):
        case=next(case for case in fixtures() if case['name']=='GMRES shifted Jordan s=4 restart 1 Jacobi 0')
        bad=copy.deepcopy(self.packed(case));n=case['a']['rows'];length=int(bad['values'][2])
        bad['values'][3+n+length+1]*=.5
        with self.assertRaisesRegex(AssertionError,'incorrect analytic residual history'):
            check_result(bad,case)
    def test_independent_pivoted_lu_on_nonsymmetric_transport(self):
        a=transport(3,1,.2,4,30);expected=[float(i+1) for i in range(9)]
        result=dense_solve(a,multiply(a,expected))
        for x,y in zip(result,expected):self.assertAlmostEqual(x,y,places=12)
        entries={(i,a['indices'][p]):a['values'][p] for i in range(9) for p in range(a['offsets'][i],a['offsets'][i+1])}
        self.assertNotEqual(entries[0,1],entries[1,0])
    def test_public_metadata_keeps_work_and_drops_secrets(self):
        row=dict(implementation='c',solver_iterations=12,restart=5,restart_count=2,logical_workspace_bytes=4096,environment={'TOKEN':'secret'},hostname='private')
        result=PublicSanitizer().report({'results':[row]})['results'][0]
        for key in ['solver_iterations','restart','restart_count','logical_workspace_bytes']:self.assertEqual(result[key],row[key])
        self.assertNotIn('environment',result);self.assertNotIn('hostname',result)

    def test_workspace_caps_restart_by_problem_and_iteration_limit(self):
        from gmres_bench import workspace_bytes
        self.assertEqual(workspace_bytes(8,100,1000),workspace_bytes(8,8,1000))
        self.assertEqual(workspace_bytes(64,20,5),workspace_bytes(64,5,1000))
        self.assertGreater(workspace_bytes(64,20,1000),workspace_bytes(64,5,1000))
