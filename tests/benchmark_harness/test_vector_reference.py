"""Check independent quaternion and exact cross-product references and workload shapes."""
import math
from pathlib import Path
import sys
from unittest import TestCase

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmarks'))
from reference import determinant, matrix
from vector_reference import vector_checksum, vector_expected, vector_fixtures, vector_inputs


class VectorReferenceTests(TestCase):
    def test_cross_oracle_uses_right_hand_rule_and_left_orientation(self):
        a = matrix(1,3,[1,2,3])
        b = matrix(3,1,[-2,4,1])
        self.assertEqual(vector_expected('cross',a,b), matrix(1,3,[-10,-7,8]))
        self.assertEqual(vector_expected('cross',b,a), matrix(3,1,[10,7,-8]))
        self.assertEqual(vector_expected('cross',a,a)['values'], [0,0,0])
        self.assertEqual(vector_expected('cross',matrix(3,1,[1,0,0]),matrix(3,1,[0,1,0]))['values'], [0,0,1])

    def test_quaternion_rotation_has_active_right_handed_orientation(self):
        expected = {
            (1,0,0): [1,0,0,0,0,-1,0,1,0],
            (0,1,0): [0,0,1,0,1,0,-1,0,0],
            (0,0,1): [0,-1,0,1,0,0,0,0,1],
        }
        for axis, wanted in expected.items():
            actual = vector_expected('rotation3d',matrix(3,1,axis),scalar=math.pi/2)['values']
            for value, reference in zip(actual,wanted):
                self.assertAlmostEqual(value,reference,places=14)
        clockwise = vector_expected('rotation2d',matrix(0,0,[]),scalar=-math.pi/2)['values']
        for value, wanted in zip(clockwise,[0,1,-1,0]):
            self.assertAlmostEqual(value,wanted,places=14)

    def test_arbitrary_axis_preserves_axis_norm_orientation_and_inverse(self):
        for axis in ([1,2,3], [1e300,2e300,3e300], [1e-300,2e-300,3e-300], [5e-324,0,0], [1.7e308]*3):
            for angle in (0,.5,-2.3,1e-10):
                a = matrix(3,1,axis)
                r = vector_expected('rotation3d',a,scalar=angle)
                inverse = vector_expected('rotation3d',a,scalar=-angle)['values']
                values = r['values']
                largest = max(map(abs,axis))
                unit = [value/largest for value in axis]
                norm = math.hypot(*unit)
                unit = [value/norm for value in unit]
                for i in range(3):
                    self.assertAlmostEqual(sum(values[i*3+k]*unit[k] for k in range(3)),unit[i],places=13)
                    for j in range(3):
                        self.assertAlmostEqual(sum(values[k*3+i]*values[k*3+j] for k in range(3)),int(i==j),places=13)
                        self.assertAlmostEqual(sum(values[i*3+k]*inverse[k*3+j] for k in range(3)),int(i==j),places=13)
                self.assertAlmostEqual(determinant(r),1,places=13)

    def test_workload_shapes_angles_and_checksum_consumption(self):
        self.assertNotEqual(vector_checksum('cross',3,17),0)
        for op,size in (('cross',3),('rotation2d',2),('rotation3d',3)):
            for seed in (0,17,2147483646):
                a,b,angle = vector_inputs(op,size,seed)
                self.assertEqual(angle,.5)
                wanted = vector_expected(op,a,b,angle)['values']
                self.assertEqual(vector_checksum(op,size,seed),wanted[0]+wanted[len(wanted)//2]+wanted[-1])
                if op == 'cross':
                    self.assertEqual((a['rows'],a['cols'],b['rows'],b['cols']),(3,1,3,1))
                    self.assertEqual(a['values'],[((i*17+seed*13)%101-50)/16 for i in range(3)])
                    self.assertEqual(b['values'][2],-((2*17+(seed+1)*13)%101-50)/16)
                elif op == 'rotation2d':
                    self.assertEqual(a,matrix(0,0,[]))
                else:
                    self.assertEqual(a,matrix(3,1,[1,2,3]))

    def test_invalid_workload_sizes_seeds_and_shapes_are_rejected(self):
        for op,size in (('cross',2),('rotation2d',3),('rotation3d',2),('unknown',3),('cross',True)):
            with self.assertRaises(ValueError): vector_inputs(op,size,0)
        for seed in (-1,2147483647,.5,True):
            with self.assertRaises(ValueError): vector_inputs('cross',3,seed)
        for case in vector_fixtures():
            if case['invalid']:
                with self.subTest(name=case['name']), self.assertRaises((ValueError,OverflowError)):
                    vector_expected(case['op'],case['a'],case['b'],case['scalar'])

    def test_fixtures_are_unique_finite_and_cover_orientations_and_scales(self):
        cases = vector_fixtures()
        self.assertEqual(len(cases),len({case['name'] for case in cases}))
        self.assertEqual({case['op'] for case in cases},{'cross','rotation2d','rotation3d'})
        self.assertGreaterEqual(sum(case['invalid'] for case in cases),10)
        for case in cases:
            for operand in (case['a'],case['b']):
                if operand is not None:
                    self.assertTrue(all(math.isfinite(value) for value in operand['values']))
            if not case['invalid']:
                result = case['expected']
                self.assertEqual(len(result['values']),result['rows']*result['cols'])
                self.assertTrue(all(math.isfinite(value) for value in result['values']))
