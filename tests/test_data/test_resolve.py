"""Tests the tp.data.resolve module

This module tests whether the resolve function resolves the correct
array depth and temperature. For tests on the direction, see
test_data/test_aniso.
Failures in tp.data.aniso can cause knock-on effects here.
Failures here can cause knock on effect in tp.plot modules.
"""

import unittest
from tp.data import resolve

class ResolveTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.one = [1., 2., 3.]
        cls.two = [[ 1.,  2.,  3.],
                   [11., 12., 13.],
                   [21., 22., 23.]]
        cls.three = [[[ 1.,  2.,  3.]],
                     [[11., 12., 13.]],
                     [[21., 22., 23.]]]
        cls.four = [[[[ 1.,  2.,  3.]]],
                    [[[11., 12., 13.]]],
                    [[[21., 22., 23.]]]]
        cls.three_m = [[[[ 1.,  2.,  3.], [ 4.,  5.,  6.], [ 7.,  8.,  9.]],
                        [[11., 12., 13.], [14., 15., 16.], [17., 18., 19.]],
                        [[21., 22., 23.], [24., 25., 26.], [27., 28., 29.]]]]
        cls.t = [10., 20., 30.]

    def test_conductivity_direction(self):
        q = 'conductivity'
        d = {q:              self.three_m,
             'temperature':  self.t,
             'meta':         {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [[[9., 19., 29.]]]).all())

    def test_conductivity_temperature(self):
        q = 'conductivity'
        d = {q:              self.three_m,
             'temperature':  self.t,
             'meta':         {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.three_m,
                        'If temperature resolution is now implemented, '
                        'this test will fail, edit it!')
        # Once temperature resolution in amset is implemented, use the
        # following block
        #self.assertEqual(d[q], [[11., 12., 13.],
        #                        [14., 15., 16.],
        #                        [17., 18., 19.]])

    def test_electronic_thermal_conductivity_direction(self):
        q = 'electronic_thermal_conductivity'
        d = {q:             self.three_m,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [[[9., 19., 29.]]]).all())

    def test_electronic_thermal_conductivity_temperature(self):
        q = 'electronic_thermal_conductivity'
        d = {q:             self.three_m,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.three_m,
                        'If temperature resolution is now implemented, '
                        'this test will fail, edit it!')
        # Once temperature resolution in amset is implemented, use the
        # following block
        #self.assertEqual(d[q], [[11., 12., 13.],
        #                        [14., 15., 16.],
        #                        [17., 18., 19.]])

    def test_gamma_direction(self):
        q = 'gamma'
        d = {q:             self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertEqual(d[q], self.three)

    def test_gamma_temperature(self):
        q = 'gamma'
        d = {q      :       self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], [[11., 12., 13.]])

    def test_group_velocity_direction(self):
        q = 'group_velocity'
        d = {q:             self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [[[3.], [13.], [23.]]]).all())

    def test_group_velocity_temperature(self):
        q = 'group_velocity'
        d = {q:             self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.three)

    def test_gv_by_gv_direction(self):
        q = 'gv_by_gv'
        d = {q:             self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [[[3.], [13.], [23.]]]).all())

    def test_gv_by_gv_temperature(self):
        q = 'gv_by_gv'
        d = {q:             self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.three)

    def test_heat_capacity_direction(self):
        q = 'heat_capacity'
        d = {q:             self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertEqual(d[q], self.three)

    def test_heat_capacity_temperature(self):
        q = 'heat_capacity'
        d = {q      :       self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], [[11., 12., 13.]])

    def test_lattice_thermal_conductivity_direction(self):
        q = 'lattice_thermal_conductivity'
        d = {q:             self.two,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [3., 13., 23.]).all())

    def test_lattice_thermal_conductivity_temperature(self):
        q = 'lattice_thermal_conductivity'
        d = {q      :       self.two,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], [11., 12., 13.])

    def test_lattice_thermal_conductivity_both(self):
        q = 'lattice_thermal_conductivity'
        d = {q      :       self.two,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23., direction='z')
        self.assertEqual(d[q], 13.)

    def test_lifeftime_direction(self):
        q = 'lifetime'
        d = {q:             self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertEqual(d[q], self.three)

    def test_lifetime_temperature(self):
        q = 'lifetime'
        d = {q      :       self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], [[11., 12., 13.]])

    def test_mean_free_path_direction(self):
        q = 'mean_free_path'
        d = {q:             self.two,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertEqual(d[q], self.two)

    def test_mean_free_path_temperature(self):
        q = 'mean_free_path'
        d = {q      :       self.two,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.two)

    def test_mesh_direction(self):
        q = 'mesh'
        d = {q:             self.one,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertEqual(d[q], 3)

    def test_mesh_temperature(self):
        q = 'mesh'
        d = {q:             self.one,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.one)

    def test_mode_kappa_direction(self):
        q = 'mode_kappa'
        d = {q:             self.four,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [[[[3.]], [[13.]], [[23.]]]]).all())

    def test_mode_kappa_temperature(self):
        q = 'mode_kappa'
        d = {q      :       self.four,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], [[[11., 12., 13.]]])

    def test_mode_kappa_both(self):
        q = 'mode_kappa'
        d = {q      :       self.four,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23., direction='z')
        self.assertEqual(d[q], [[[13.]]])

    def test_occupation_direction(self):
        q = 'occupation'
        d = {q:             self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertEqual(d[q], self.three)

    def test_occupation_temperature(self):
        q = 'occupation'
        d = {q      :       self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], [[11., 12., 13.]])

    def test_power_factor_direction(self):
        q = 'power_factor'
        d = {q:             self.three_m,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [[[9., 19., 29.]]]).all())

    def test_power_factor_temperature(self):
        q = 'power_factor'
        d = {q:             self.three_m,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.three_m,
                        'If temperature resolution is now implemented, '
                        'this test will fail, edit it!')
        # Once temperature resolution in amset is implemented, use the
        # following block
        #self.assertEqual(d[q], [[11., 12., 13.],
        #                        [14., 15., 16.],
        #                        [17., 18., 19.]])

    def test_qpoint_direction(self):
        q = 'qpoint'
        d = {q:             self.two,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [3., 13., 23.]).all())

    def test_qpoint_temperature(self):
        q = 'qpoint'
        d = {q:             self.three,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.three)

    def test_seebeck_direction(self):
        q = 'seebeck'
        d = {q:             self.three_m,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [[[9., 19., 29.]]]).all())

    def test_seebeck_temperature(self):
        q = 'seebeck'
        d = {q:             self.three_m,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.three_m,
                        'If temperature resolution is now implemented, '
                        'this test will fail, edit it!')
        # Once temperature resolution in amset is implemented, use the
        # following block
        #self.assertEqual(d[q], [[11., 12., 13.],
        #                        [14., 15., 16.],
        #                        [17., 18., 19.]])

    def test_zt_direction(self):
        q = 'zt'
        d = {q:             self.three_m,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, direction='z')
        self.assertTrue((d[q] == [[[9., 19., 29.]]]).all())

    def test_zt_temperature(self):
        q = 'zt'
        d = {q:             self.three_m,
             'temperature': self.t,
             'meta':        {}}
        d = resolve.resolve(d, q, temperature=23.)
        self.assertEqual(d[q], self.three_m,
                        'If temperature resolution is now implemented, '
                        'this test will fail, edit it!')
        # Once temperature resolution in amset is implemented, use the
        # following block
        #self.assertEqual(d[q], [[11., 12., 13.],
        #                        [14., 15., 16.],
        #                        [17., 18., 19.]])

if __name__ == '__main__':
    unittest.main()
