import json
import math
from pathlib import Path
import subprocess
import sys
import unittest

from matrix import Matrix


class VectorRotationTests(unittest.TestCase):
    def assert_matrix_close(self, actual, expected, tolerance=1e-13):
        self.assertEqual((actual.rows,actual.cols),(expected.rows,expected.cols))
        for value,wanted in zip(actual.values,expected.values):
            self.assertLessEqual(abs(value-wanted),tolerance)

    def test_cross_is_right_handed_preserves_left_shape_and_does_not_mutate(self):
        for left_shape in ((3,1),(1,3)):
            for right_shape in ((3,1),(1,3)):
                a = Matrix(*left_shape,[1,2,3])
                b = Matrix(*right_shape,[-2,4,1])
                result = a.cross(b)
                self.assertEqual((result.rows,result.cols),left_shape)
                self.assertEqual(result.values,[-10,-7,8])
                self.assertEqual(b.cross(a).values,[10,7,-8])
                self.assertEqual(a.cross(a).values,[0,0,0])
                result[0,0] = 99
                self.assertEqual(a.values,[1,2,3])
                self.assertEqual(b.values,[-2,4,1])
        self.assertEqual(Matrix(3,1,[1,0,0]).cross(Matrix(3,1,[0,1,0])).values,[0,0,1])
        self.assertEqual(Matrix(3,1,[1e300,0,0]).cross(Matrix(1,3,[0,1e-300,0])).values,[0,0,1])

    def test_cross_is_orthogonal_bilinear_and_satisfies_lagrange_identity(self):
        a,b,c = Matrix(3,1,[.5,-1.25,2.5]),Matrix(3,1,[3,.25,-2]),Matrix(3,1,[2,-.5,.25])
        cross = a.cross(b).values
        dot = lambda left,right: sum(x*y for x,y in zip(left,right))
        self.assertAlmostEqual(dot(cross,a.values),0)
        self.assertAlmostEqual(dot(cross,b.values),0)
        self.assertAlmostEqual(dot(cross,cross),dot(a.values,a.values)*dot(b.values,b.values)-dot(a.values,b.values)**2)
        self.assert_matrix_close(a.add(c).cross(b),a.cross(b).add(c.cross(b)))
        self.assert_matrix_close(a.scale(-.5).cross(b),a.cross(b).scale(-.5))

    def test_principal_rotations_have_the_column_vector_right_hand_convention(self):
        self.assert_matrix_close(Matrix.rotation2d(math.pi/2).multiply(Matrix(2,1,[1,0])),Matrix(2,1,[0,1]))
        self.assert_matrix_close(Matrix.rotation2d(-math.pi/2).multiply(Matrix(2,1,[1,0])),Matrix(2,1,[0,-1]))
        cases = ((Matrix.rotation_x,[0,1,0],[0,0,1],[1,0,0]),
                 (Matrix.rotation_y,[0,0,1],[1,0,0],[0,1,0]),
                 (Matrix.rotation_z,[1,0,0],[0,1,0],[0,0,1]))
        for rotation,basis,wanted,axis in cases:
            result = rotation(math.pi/2)
            self.assert_matrix_close(result.multiply(Matrix(3,1,basis)),Matrix(3,1,wanted))
            self.assert_matrix_close(result,Matrix.rotation_axis_angle(Matrix(1,3,axis),math.pi/2))

    def test_arbitrary_axis_rotations_are_proper_orthogonal_and_compose(self):
        for entries in ([1,2,3],[1e300,2e300,3e300],[1e-300,2e-300,3e-300],[5e-324,0,0],[1.7e308]*3):
            axis = Matrix(3,1,entries)
            before = axis.values
            r = Matrix.rotation_axis_angle(axis,.5)
            self.assertEqual(axis.values,before)
            self.assert_matrix_close(r.transpose().multiply(r),Matrix.identity(3))
            self.assertAlmostEqual(r.determinant(),1,places=13)
            self.assert_matrix_close(r.transpose(),Matrix.rotation_axis_angle(axis,-.5))
            self.assert_matrix_close(r.transpose(),Matrix.rotation_axis_angle(axis.scale(-1),.5))
            self.assert_matrix_close(r.multiply(Matrix.rotation_axis_angle(axis,.7)),Matrix.rotation_axis_angle(axis,1.2))
            largest = max(map(abs,entries))
            direction = Matrix(3,1,[value/largest for value in entries])
            self.assert_matrix_close(r.multiply(direction),direction)
        self.assert_matrix_close(Matrix.rotation_axis_angle(Matrix(3,1,[1,2,3]),0),Matrix.identity(3))
        tiny = Matrix.rotation_axis_angle(Matrix(3,1,[1,1,0]),1e-10)
        self.assertGreater(tiny[0,1],0)
        self.assertAlmostEqual(tiny[0,1]/2.5e-21,1,places=13)

    def test_invalid_vectors_angles_axes_and_cross_overflow_fail_without_mutation(self):
        a = Matrix(3,1,[1,2,3])
        for wrong in (Matrix(),Matrix(2,1,[1,2]),Matrix(2,2,[1,2,3,4])):
            with self.assertRaises(ValueError): a.cross(wrong)
            with self.assertRaises(ValueError): wrong.cross(a)
            with self.assertRaises(ValueError): Matrix.rotation_axis_angle(wrong,.5)
        with self.assertRaises(ValueError): Matrix.rotation_axis_angle(Matrix(3,1),.5)
        with self.assertRaises(ValueError): a.cross(None)
        for angle in (math.nan,math.inf,-math.inf):
            for rotation in (Matrix.rotation2d,Matrix.rotation_x,Matrix.rotation_y,Matrix.rotation_z):
                with self.assertRaises(ValueError): rotation(angle)
            with self.assertRaises(ValueError): Matrix.rotation_axis_angle(a,angle)
        large = Matrix(3,1,[1e308,0,0])
        with self.assertRaises(ValueError): large.cross(Matrix(3,1,[0,1e308,0]))
        self.assertEqual(large.values,[1e308,0,0])
        self.assertEqual(a.values,[1,2,3])

    def test_vector_rotation_runner_protocol_and_fixed_benchmark_sizes(self):
        runner = str(Path(__file__).with_name('runner.py'))
        def invoke(args,data=''):
            return subprocess.run([sys.executable,runner]+args,input=data,text=True,capture_output=True)
        result = invoke(['check','cross','1','3','3','1'],'1 2 3 -2 4 1')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(json.loads(result.stdout),{'rows':1,'cols':3,'values':[-10,-7,8]})
        for op,args,data,dimension in (('rotation2d',['0','0',str(math.pi/2)],'',2),
                                      ('rotation3d',['3','1',str(math.pi/2)],'0 0 1',3)):
            result = invoke(['check',op]+args,data)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(json.loads(result.stdout)['rows'],dimension)
        for op,size in (('cross',3),('rotation2d',2),('rotation3d',3)):
            result = invoke(['bench',op,str(size),'2','17'])
            self.assertEqual(result.returncode,0,result.stderr)
            payload = json.loads(result.stdout)
            self.assertTrue(math.isfinite(payload['checksum']))
            self.assertNotEqual(payload['checksum'],0)
            self.assertGreater(payload['elapsed_ns'],0)
            self.assertNotEqual(invoke(['bench',op,str(size+1),'2','17']).returncode,0)
        self.assertNotEqual(invoke(['check','rotation2d','1','1','.5'],'1').returncode,0)
        self.assertNotEqual(invoke(['check','rotation2d','0','0','.5'],'1').returncode,0)
