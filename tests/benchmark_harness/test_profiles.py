"""Profile parsers preserve observed stacks and avoid inclusive-count double counting."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks"))
from profiles import parse_go_raw, parse_macos_sample, parse_node, parse_python


class ProfileParserTests(unittest.TestCase):
    def test_node_resolves_parents_and_microsecond_deltas(self):
        stacks, weights = parse_node({"nodes": [
            {"id": 1, "callFrame": {"functionName": "(root)"}, "children": [2]},
            {"id": 2, "callFrame": {"functionName": "multiply", "url": "file:///matrix.js", "lineNumber": 20}},
        ], "samples": [2, 1], "timeDeltas": [1500, 500]})
        self.assertEqual(stacks, [["(root)", "multiply (matrix.js:21)"], ["(root)"]])
        self.assertEqual(weights, [1.5, 0.5])

    def test_python_uses_real_sample_intervals(self):
        stacks, weights = parse_python({"samples": [
            {"at_ns": 1000000, "stack": ["main", "multiply"]},
            {"at_ns": 7500000, "stack": ["main", "cleanup"]},
        ], "end_ns": 9000000})
        self.assertEqual(stacks, [["main", "multiply"], ["main", "cleanup"]])
        self.assertEqual(weights, [6.5, 1.5])

    def test_macos_converts_inclusive_tree_to_exclusive_stack_weights(self):
        text = """Analysis of sampling runner (pid 123) every 1 millisecond
Call graph:
    10 Thread_123 DispatchQueue_1: com.apple.main-thread (serial)
      10 main (in runner) + 32 [0x10]
        7 multiply (in runner) + 20 [0x20]
        + 5 dot (in runner) + 12 [0x30]
        2 allocate (in runner) + 4 [0x40]
    10 Thread_999 worker
      10 unrelated (in system) + 0 [0x50]
Total number in stack:
"""
        stacks, weights = parse_macos_sample(text)
        self.assertEqual(sum(weights), 10)
        self.assertEqual(weights, [1, 2, 5, 2])
        self.assertEqual(stacks[2], ["main thread", "main (in runner)", "multiply (in runner)", "dot (in runner)"])
        self.assertFalse(any("unrelated" in frame for stack in stacks for frame in stack))

    def test_go_reverses_leaf_first_stacks_and_keeps_inline_frames(self):
        text = """PeriodType: cpu nanoseconds
Samples:
samples/count cpu/nanoseconds
  3 30000000: 1 2
Locations
  1: 0x100 M=1 math.IsInf /src/math.go:1:0 s=1
          matrix.Multiply /src/matrix.go:10:0 s=2
  2: 0x200 M=1 main.main /src/main.go:20:0 s=3
Mappings
1: ignored
"""
        stacks, weights = parse_go_raw(text)
        self.assertEqual(weights, [30])
        self.assertTrue(stacks[0][0].startswith("main.main"))
        self.assertTrue(stacks[0][1].startswith("matrix.Multiply"))
        self.assertTrue(stacks[0][2].startswith("math.IsInf"))


if __name__ == "__main__":
    unittest.main()
