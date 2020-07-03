"""Tests the tp.data.aniso module

Failures here can have knock-on effects in tp.data.resolve and tp.plot
modules.
"""

import unittest
from tp.data import aniso

class OneTest(unittest.TestCase):
    def setUp(self):
        self.a = [2., 5., 14.]
        self.func = aniso.one

    def test_x(self):
        self.assertEqual(self.func(self.a, 'x'), 2.)

    def test_y(self):
        self.assertEqual(self.func(self.a, 'y'), 5.)

    def test_z(self):
        self.assertEqual(self.func(self.a, 'z'), 14.)

    def test_average(self):
        self.assertEqual(self.func(self.a, 'average'), 7)

    def test_norm(self):
        self.assertEqual(self.func(self.a, 'norm'), 15)

class TwoTest(unittest.TestCase):
    def setUp(self):
        self.a = [[2., 5., 14.]]
        self.func = aniso.two

    def test_x(self):
        self.assertEqual(self.func(self.a, 'x'), [2.])

    def test_y(self):
        self.assertEqual(self.func(self.a, 'y'), [5.])

    def test_z(self):
        self.assertEqual(self.func(self.a, 'z'), [14.])

    def test_average(self):
        self.assertEqual(self.func(self.a, 'average'), [7])

    def test_norm(self):
        self.assertEqual(self.func(self.a, 'norm'), [15])

class ThreeTest(unittest.TestCase):
    def setUp(self):
        self.a = [[[2., 5., 14.]]]
        self.func = aniso.three

    def test_x(self):
        self.assertEqual(self.func(self.a, 'x'), [[2.]])

    def test_y(self):
        self.assertEqual(self.func(self.a, 'y'), [[5.]])

    def test_z(self):
        self.assertEqual(self.func(self.a, 'z'), [[14.]])

    def test_average(self):
        self.assertEqual(self.func(self.a, 'average'), [[7]])

    def test_norm(self):
        self.assertEqual(self.func(self.a, 'norm'), [[15]])

class FourTest(unittest.TestCase):
    def setUp(self):
        self.a = [[[[2., 5., 14.]]]]
        self.func = aniso.four

    def test_x(self):
        self.assertEqual(self.func(self.a, 'x'), [[[2.]]])

    def test_y(self):
        self.assertEqual(self.func(self.a, 'y'), [[[5.]]])

    def test_z(self):
        self.assertEqual(self.func(self.a, 'z'), [[[14.]]])

    def test_average(self):
        self.assertEqual(self.func(self.a, 'average'), [[[7]]])

    def test_norm(self):
        self.assertEqual(self.func(self.a, 'norm'), [[[15]]])

class MatrixThreeTest(unittest.TestCase):
    def setUp(self):
        self.a = [[[[2., 3., 4.],
                    [1., 5., 6.],
                    [8., 9., 14.]]]]
        self.func = aniso.matrix_three

    def test_x(self):
        self.assertEqual(self.func(self.a, 'x'), [[2.]])

    def test_y(self):
        self.assertEqual(self.func(self.a, 'y'), [[5.]])

    def test_z(self):
        self.assertEqual(self.func(self.a, 'z'), [[14.]])

    def test_average(self):
        self.assertEqual(self.func(self.a, 'average'), [[7]])

    def test_norm(self):
        self.assertEqual(self.func(self.a, 'norm'), [[15]])

if __name__ == '__main__':
    unittest.main()
