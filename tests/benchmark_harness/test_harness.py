"""Tests for the independent oracle, checksum gate, and offline report encoding."""
import json
import math
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"benchmarks"))
import reference
import report
import run


class ReferenceTests(unittest.TestCase):
    def test_exact_oracle_handles_sign_empty_and_fractional_values(self):
        self.assertEqual(reference.determinant(reference.matrix(0,0,[])),1)
        self.assertEqual(reference.determinant(reference.matrix(3,3,[6,1,1,4,-2,5,2,8,7])),-306)
        self.assertEqual(reference.determinant(reference.matrix(2,2,[0.5,1.25,-2,0.75])),2.875)
        self.assertEqual(reference.determinant(reference.matrix(2,2,[3,4,1,2])),2)

    def test_benchmark_checksum_matches_full_independent_results(self):
        for size in (1,2,4):
            for op in reference.OPERATIONS:
                a=reference.generated(size,17,op=="determinant")
                b=reference.generated(size,18)
                result=reference.expected(op,a,b)
                if "value" in result:
                    checksum=result["value"]
                else:
                    values=result["values"]
                    checksum=values[0]+values[len(values)//2]+values[-1]
                self.assertTrue(math.isclose(reference.benchmark_checksum(op,size,17),float(checksum),rel_tol=1e-12))

    def test_wrong_shape_and_nonfinite_outputs_are_rejected(self):
        for bad in ({"rows":3,"cols":2,"values":[1]*6}, {"rows":2,"cols":3,"values":[math.nan]*6}):
            with self.assertRaises(AssertionError):
                reference.assert_result(bad,reference.matrix(2,3,[1]*6))

    def test_fixture_names_are_unique(self):
        cases=reference.fixtures()
        self.assertEqual(len(cases),len({c["name"] for c in cases}))
        self.assertGreaterEqual(sum(c["invalid"] for c in cases),5)

    def test_malformed_result_types_never_pass(self):
        expected=reference.matrix(1,1,[1])
        for bad in (None,7,[],{"rows":True,"cols":1,"values":[1]},
                    {"rows":1.0,"cols":1,"values":[1]},
                    {"rows":1,"cols":1,"values":None},
                    {"rows":1,"cols":1,"values":"1"},
                    {"rows":1,"cols":1,"values":[True]},
                    {"rows":1,"cols":1,"values":[10**400]}):
            with self.subTest(result=bad), self.assertRaises(AssertionError):
                reference.assert_result(bad,expected)
        with self.assertRaises(AssertionError):
            reference.assert_result({"upper":1,"lower":0},{"upper":True,"lower":False})
        with self.assertRaises(AssertionError):
            reference.assert_result({"value":False},{"value":0})


class MeasurementTests(unittest.TestCase):
    def setUp(self):
        self.implementation={"runner":["unused"],"env":{}}

    def test_corrupted_checksum_never_becomes_a_timing(self):
        with patch.object(run,"require",return_value=json.dumps({"elapsed_ns":200,"iterations":2,"checksum":999})):
            with self.assertRaises(ValueError): run.measure(self.implementation,"add",2,2,17,3)

    def test_zero_timing_wrong_iterations_and_nonfinite_checksum_rejected(self):
        for payload in ({"elapsed_ns":0,"iterations":2,"checksum":6},
                        {"elapsed_ns":100,"iterations":3,"checksum":6},
                        {"elapsed_ns":100,"iterations":2,"checksum":float("inf")}):
            with patch.object(run,"require",return_value=json.dumps(payload)):
                with self.assertRaises(ValueError): run.measure(self.implementation,"add",2,2,17,3)

    def test_per_operation_timing_uses_reported_elapsed_and_known_count(self):
        with patch.object(run,"require",return_value='{"elapsed_ns":100,"iterations":2,"checksum":6}'):
            self.assertEqual(run.measure(self.implementation,"add",2,2,17,3)["ns_per_op"],50)

    def test_eigen_benchmark_rejects_wrong_basis_before_starting_clock(self):
        # An identity basis has the right norm and therefore the same checksum,
        # but cannot diagonalize these tridiagonal inputs.
        def wrong_basis(command, **kwargs):
            self.assertEqual(command[1:3], ["check", "eigen_symmetric"])
            size = int(command[3])
            self.assertEqual(len(kwargs["stdin"].split()), size * size)
            return json.dumps({"eigenvalues": run.eigen_spectrum(size, 17),
                               "eigenvectors": reference.matrix(size, size,
                                   [float(row == col) for row in range(size) for col in range(size)])})

        implementation = dict(self.implementation, id="test")
        options = SimpleNamespace(seed=17, suite="quick", samples=3, operations=["eigen_symmetric"])
        with patch.object(run, "require", side_effect=wrong_basis) as probe, \
             patch.object(run, "measure") as measure, patch("builtins.print"):
            results = run.benchmark([implementation], options)
        self.assertEqual(probe.call_count, 2)
        measure.assert_not_called()
        self.assertEqual([row["size"] for row in results], [8, 16])
        self.assertTrue(all(row["status"] == "failed" and row["samples"] == [] for row in results))

    def test_vector_benchmarks_check_full_results_before_starting_clock(self):
        for op, size in (("cross",3), ("rotation2d",2), ("rotation3d",3)):
            a,b,scalar = run.vector_inputs(op,size,17)
            wrong = run.vector_expected(op,a,b,scalar)
            # Corrupt two consumed entries without changing their checksum sum.
            wrong["values"][0] += 1
            wrong["values"][-1] -= 1
            options = SimpleNamespace(seed=17,suite="quick",samples=3,operations=[op])
            with self.subTest(operation=op), patch.object(run,"require",return_value=json.dumps(wrong)) as probe, \
                 patch.object(run,"measure") as measure, patch("builtins.print"):
                results = run.benchmark([dict(self.implementation,id="test")],options)
            self.assertEqual(probe.call_count,1)
            measure.assert_not_called()
            self.assertEqual(results[0]["status"],"failed")

    def test_malformed_timing_types_never_pass(self):
        for payload in (None,[],"timing",{"elapsed_ns":100,"iterations":1},
                        {"elapsed_ns":100,"iterations":True,"checksum":0},
                        {"elapsed_ns":100,"iterations":1.0,"checksum":0},
                        {"elapsed_ns":100,"iterations":1,"checksum":False},
                        {"elapsed_ns":True,"iterations":1,"checksum":0},
                        {"elapsed_ns":10**400,"iterations":1,"checksum":0},
                        {"elapsed_ns":100,"iterations":1,"checksum":10**400},
                        {"elapsed_ns":100,"iterations":1,"checksum":0,"extra":0}):
            with self.subTest(payload=payload), patch.object(run,"require",return_value=json.dumps(payload)):
                with self.assertRaises(ValueError):
                    run.measure(self.implementation,"trace",1,1,17,0)


class VerificationTests(unittest.TestCase):
    def test_malformed_output_is_a_failed_check(self):
        case={"name":"scalar","op":"trace","a":reference.matrix(1,1,[1]),
              "b":None,"invalid":False,"expected":{"value":1}}
        for output in ("null",'{"value":false}'):
            implementation={"runner":["unused"],"env":{},"status":"built"}
            completed=SimpleNamespace(returncode=0,stdout=output,stderr="")
            with patch.object(run,"fixtures",return_value=[case]), patch.object(run,"execute",return_value=completed):
                self.assertFalse(run.verify(implementation))
            self.assertEqual(implementation["status"],"failed")
            self.assertEqual(len(implementation["checks"]),1)
            self.assertFalse(implementation["checks"][0]["passed"])

    def test_missing_built_runner_is_a_failed_check(self):
        implementation={"runner":["missing"],"env":{},"status":"built"}
        with patch.object(run,"fixtures",return_value=reference.fixtures()[:1]), \
             patch.object(run,"execute",side_effect=FileNotFoundError("missing runner")):
            self.assertFalse(run.verify(implementation))
        self.assertEqual(implementation["status"],"failed")

    def test_invalid_request_must_keep_stdout_empty(self):
        invalid=next(case for case in reference.fixtures() if case["invalid"])
        completed=SimpleNamespace(returncode=1,stdout='{"value":0}',stderr="bad shape")
        implementation={"runner":["unused"],"env":{}}
        with patch.object(run,"fixtures",return_value=[invalid]), patch.object(run,"execute",return_value=completed):
            self.assertFalse(run.verify(implementation))


class StatusTests(unittest.TestCase):
    def test_each_requested_implementation_gets_one_final_record(self):
        cases=((run.ToolchainUnavailable("missing compiler"),None,"unavailable"),
               (RuntimeError("compiler failed"),None,"failed"),
               (FileNotFoundError("missing build input"),None,"failed"),
               (None,TypeError("malformed response"),"failed"))
        for build_error,verify_error,status in cases:
            with self.subTest(status=status,error=build_error or verify_error), tempfile.TemporaryDirectory() as directory:
                implementation={"id":"python","name":"Python","status":"built","runner":["unused"],"env":{}}
                argv=["run.py","--implementations","python","--verify-only","--output",directory]
                with patch.object(sys,"argv",argv), \
                     patch.object(run,"build_one",return_value=implementation,side_effect=build_error), \
                     patch.object(run,"verify",side_effect=verify_error), \
                     patch.object(run,"require",side_effect=["revision\n",""]), \
                     patch.object(run,"fingerprint",return_value="test-fingerprint"), \
                     patch.object(report,"render"), patch("builtins.print"):
                    self.assertEqual(run.main(),1)
                data=json.loads((Path(directory)/"results.json").read_text())
                self.assertEqual(len(data["implementations"]),1)
                self.assertEqual(data["implementations"][0]["id"],"python")
                self.assertEqual(data["implementations"][0]["status"],status)

    def test_invalid_compiler_override_is_unavailable(self):
        with patch.dict(run.os.environ,{"GO":"/nonexistent/linear-a-go"}), \
             patch.object(run.shutil,"which",return_value=None):
            with self.assertRaises(run.ToolchainUnavailable):
                run.binary("GO","go")


class ReportTests(unittest.TestCase):
    def test_profile_symbols_cannot_end_the_embedded_data_script(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"report.html"
            report.render({"profiles":[{"stacks":[["</script><script>alert('x')</script>"]],
                                          "weights":[1]}]},path)
            html=path.read_text()
            self.assertNotIn("</script><script>alert",html)
            self.assertIn('\\u003c/script>',html)
            self.assertNotIn("__REPORT_DATA__",html)


if __name__=="__main__": unittest.main()
