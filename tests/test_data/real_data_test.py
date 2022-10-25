"""Tests the tp.data.load module.

Note that the data files used originally are too large for github, but
these tests should work for any file of the correct format. It is not
automatically run by unittest, as it requires outside input.
"""

import unittest
from tp.data import load
import numpy as np
import os

class AmsetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        f = 'data/amset_data_85x85x47.json'
        cls.qs = ['doping', 'conductivity', 'electronic_thermal_conductivity',
                  'seebeck', 'temperature']
        cls.d = load.amset(f, cls.qs)
        cls.ts = len(cls.d['temperature'])
        cls.ds = len(cls.d['doping'])

    def test_conductivity(self):
        self.assertEqual(np.shape(self.d['conductivity']),
                         (self.ts, self.ds, 3, 3))

    def test_electronic_thermal_conductivity(self):
        self.assertEqual(np.shape(self.d['electronic_thermal_conductivity']),
                         (self.ts, self.ds, 3, 3))

    def test_seebeck(self):
        self.assertEqual(np.shape(self.d['seebeck']),
                         (self.ts, self.ds, 3, 3))

class BoltztrapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        f = 'data/boltztrap.hdf5'
        cls.qs = ['doping', 'conductivity', 'electronic_thermal_conductivity',
                  'seebeck', 'temperature', 'power_factor']
        cls.d = load.boltztrap(f, cls.qs)
        cls.ts = len(cls.d['temperature'])
        cls.ds = len(cls.d['doping'])

    def test_conductivity(self):
        self.assertEqual(np.shape(self.d['conductivity']),
                         (self.ts, self.ds, 3, 3))

    def test_electronic_thermal_conductivity(self):
        self.assertEqual(np.shape(self.d['electronic_thermal_conductivity']),
                         (self.ts, self.ds, 3, 3))

    def test_seebeck(self):
        self.assertEqual(np.shape(self.d['seebeck']),
                         (self.ts, self.ds, 3, 3))

    def test_power_factor(self):
        self.assertEqual(np.shape(self.d['power_factor']),
                         (self.ts, self.ds, 3, 3))

class Phono3pyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        f = 'data/kappa-m505028.hdf5'
        cls.qs = ['frequency', 'gamma', 'group_velocity',
                  'lattice_thermal_conductivity', 'lifetime', 'mean_free_path',
                  'mode_kappa', 'occupation', 'qpoint', 'temperature', 'weight']
        cls.d = load.phono3py(f, cls.qs)
        cls.ts = len(cls.d['temperature'])
        cls.qpts = len(cls.d['qpoint'])
        cls.bands = len(cls.d['frequency'][0])

    def test_frequency(self):
        self.assertEqual(np.shape(self.d['frequency']),
                         (self.qpts, self.bands))

    def test_gamma(self):
        self.assertEqual(np.shape(self.d['gamma']),
                         (self.ts, self.qpts, self.bands))

    def test_group_velocity(self):
        self.assertEqual(np.shape(self.d['group_velocity']),
                         (self.qpts, self.bands, 3))

    def test_lattice_thermal_conductivity(self):
        self.assertEqual(np.shape(self.d['lattice_thermal_conductivity']),
                         (self.ts, 6))

    def test_lifetime(self):
        self.assertEqual(np.shape(self.d['lifetime']),
                         (self.ts, self.qpts, self.bands))

    def test_mean_free_path(self):
        self.assertEqual(np.shape(self.d['mean_free_path']),
                         (self.ts, self.qpts, self.bands, 3))

    def test_mode_kappa(self):
        self.assertEqual(np.shape(self.d['mode_kappa']),
                         (self.ts, self.qpts, self.bands, 6))

    def test_occupation(self):
        self.assertEqual(np.shape(self.d['occupation']),
                         (self.ts, self.qpts, self.bands))

    def test_weight(self):
        self.assertEqual(np.shape(self.d['weight']),
                         (self.qpts,))

class PhonopyDispersionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        f = 'data/band.yaml'
        cls.d = load.phonopy_dispersion(f)
        dx = {'tick_position': list(range(len(cls.d['tick_position'])))}
        cls.d2 = load.phonopy_dispersion(f, dx)

    def test_qpoint(self):
        self.assertEqual(len(self.d['x']), len(self.d['qpoint']))

    def test_frequency(self):
        self.assertEqual(len(self.d['x']), len(self.d['frequency']))

    def test_ticks(self):
        self.assertEqual(len(self.d['tick_position']), len(self.d['tick_label']))

    def test_scale_xmin(self):
        self.assertEqual(self.d2['x'][0], 0)

    def test_scale_xmax(self):
        self.assertEqual(self.d2['x'][-1], len(self.d['tick_label']) - 1)

    def test_scale_ticks(self):
        self.assertEqual(self.d2['tick_position'],
                         list(range(len(self.d['tick_label']))))

class PhonopyDosTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        f = 'data/projected_dos.dat'
        atoms = 'Sb 2 Mg 3'
        cls.d = load.phonopy_dos(f, atoms)

    def test_dos(self):
        for d in self.d:
            if d != 'meta':
                self.assertEqual(len(self.d['frequency']), len(self.d[d]),
                                 '{} failed'.format(d))

if __name__ == '__main__':
    unittest.main()