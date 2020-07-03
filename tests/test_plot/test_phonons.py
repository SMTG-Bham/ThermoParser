"""Tests the tp.plot.phonons module

Failures in tp.plot.axes can have knock-on effects here.
"""

import unittest
from tp.plot import phonons
from tp.plot.axes import one_small_legend
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import warnings

class DispersionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'x':             [0, 1],
                 'qpoint':        [[0, 1, 2],
                                   [3, 1, 2]],
                 'eigenvalue':    [[0, 1, 2],
                                   [2, 1, 0]],
                 'tick_position': [0, 1],
                 'tick_label':    ['A', 'B'],
                 'meta':          {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ax = phonons.add_dispersion(cls.ax, cls.d, label='test',
                     colour=['red', '#00ff00'], linestyle=['-', ':'])

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (0, 1))

    def test_legend(self):
        labels = self.ax.get_legend_handles_labels()[1]
        self.assertEqual(labels, ['test'])

class AltDispersionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # warnings ignored as there are tonnes, likely due to random
        # input data to pymatgen (but filtering pymatgen warnings failed).
        warnings.simplefilter('ignore')
        poscar = 'data/POSCAR'
        cls.d = {'qpoint':      [[0, 1, 2.1],
                                 [3.2, 0.9, 2],
                                 [0.1, 0, 0.1],
                                 [1, 1, 1],
                                 [4, 1, 0]],
                 'mode_kappa':  [[[[1]], [[2]],
                                  [[2]], [[1]],
                                  [[1]], [[2]],
                                  [[2]], [[1]],
                                  [[0]], [[3]]]],
                 'temperature': [11],
                 'meta':        {}}
        cls.pd = {'x':             [0, 1, 2, 3],
                 'qpoint':        [[0, 1, 2],
                                   [3, 1, 2],
                                   [0, 0, 0],
                                   [1, 1, 1]],
                 'eigenvalue':    [[0, 1, 2],
                                   [2, 1, 0],
                                   [3, 2, 1],
                                   [1, 2, 3]],
                 'tick_position': [0, 3],
                 'tick_label':    ['A', 'B'],
                 'meta':          {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ax = phonons.add_alt_dispersion(cls.ax, cls.d, cls.pd, 'mode_kappa',
                     colour=['red', '#00ff00'], linestyle=['-', ':'],
                     interpolate=1, direction='x', poscar=poscar)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (0, 3))

    def test_qpoints(self):
        self.assertGreater(self.ax.get_ylim()[0], 0)
        self.assertLess(self.ax.get_ylim()[1], 3)

class ProjectedDispersionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # warnings ignored as there are tonnes, likely due to random
        # input data and real POSCAR.
        warnings.simplefilter('ignore')
        cls.poscar = 'data/POSCAR'
        cls.d = {'qpoint':      [[0, 1, 2.1],
                                 [3.2, 0.9, 2],
                                 [0.1, 0, 0.1],
                                 [1, 1, 1],
                                 [4, 1, 0]],
                 'mode_kappa':  [[[[1]], [[2]],
                                  [[2]], [[1]],
                                  [[1]], [[2]],
                                  [[2]], [[1]],
                                  [[0]], [[3]]]],
                 'temperature': [11],
                 'meta':        {}}
        cls.pd = {'x':             [0, 1, 2, 3],
                 'qpoint':        [[0, 1, 2],
                                   [3, 1, 2],
                                   [0, 0, 0],
                                   [1, 1, 1]],
                 'eigenvalue':    [[0, 1, 2],
                                   [2, 1, 0],
                                   [3, 2, 1],
                                   [1, 2, 3]],
                 'tick_position': [0, 3],
                 'tick_label':    ['A', 'B'],
                 'meta':          {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ax, cls.cbar = phonons.add_projected_dispersion(cls.ax, cls.d,
                              cls.pd, 'mode_kappa', interpolate=1,
                              direction='x', poscar=cls.poscar)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (0, 3))

    def test_qpoints(self):
        self.assertEqual(np.round(self.cbar.get_clim()[0], 0), 1)
        self.assertEqual(np.round(self.cbar.get_clim()[1], 0), 2)

class AltProjectedDispersionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # warnings ignored as there are tonnes, likely due to random
        # input data and real POSCAR. Could do with a P1 POSCAR, but materials.
        warnings.simplefilter('ignore')
        cls.poscar = 'data/POSCAR'
        cls.d = {'qpoint':      [[0, 1, 2.1],
                                 [3.2, 0.9, 2],
                                 [0.1, 0, 0.1],
                                 [1, 1, 1],
                                 [4, 1, 0]],
                 'mode_kappa':  [[[[1]], [[2]],
                                  [[2]], [[1]],
                                  [[1]], [[2]],
                                  [[2]], [[1]],
                                  [[0]], [[3]]]],
                 'temperature': [11],
                 'meta':        {}}
        cls.pd = {'x':             [0, 1, 2, 3],
                 'qpoint':        [[0, 1, 2],
                                   [3, 1, 2],
                                   [0, 0, 0],
                                   [1, 1, 1]],
                 'eigenvalue':    [[0, 1, 2],
                                   [2, 1, 0],
                                   [3, 2, 1],
                                   [1, 2, 3]],
                 'tick_position': [0, 3],
                 'tick_label':    ['A', 'B'],
                 'meta':          {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ax, cls.cbar = phonons.add_alt_projected_dispersion(cls.ax, cls.d,
                              cls.pd, 'mode_kappa', 'mode_kappa',
                              interpolate=1, direction='x', poscar=cls.poscar)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (0, 3))

    def test_qpoints(self):
        self.assertEqual(np.round(self.ax.get_ylim()[0], 0), 1)
        self.assertEqual(np.round(self.ax.get_ylim()[1], 0), 2)
        self.assertEqual(np.round(self.cbar.get_clim()[0], 0), 1)
        self.assertEqual(np.round(self.cbar.get_clim()[1], 0), 2)

class WidebandTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # warnings ignored as there are tonnes, likely due to random
        # input data and real POSCAR.
        warnings.simplefilter('ignore')
        cls.poscar = 'data/POSCAR'
        cls.d = {'qpoint':      [[0, 1, 2.1],
                                 [3.2, 0.9, 2],
                                 [0.1, 0, 0.1],
                                 [1, 1, 1],
                                 [4, 1, 0]],
                 'gamma':       [[[1], [2],
                                  [2], [1],
                                  [1], [2],
                                  [2], [1],
                                  [0], [3]]],
                 'temperature': [11],
                 'meta':        {}}
        cls.pd = {'x':             [0, 1, 2, 3],
                 'qpoint':        [[0, 1, 2],
                                   [3, 1, 2],
                                   [0, 0, 0],
                                   [1, 1, 1]],
                 'eigenvalue':    [[0, 1, 2],
                                   [2, 1, 0],
                                   [3, 2, 1],
                                   [1, 2, 3]],
                 'tick_position': [0, 3],
                 'tick_label':    ['A', 'B'],
                 'meta':          {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ax = phonons.add_wideband(cls.ax, cls.d, cls.pd, interpolate=1,
                     poscar=cls.poscar)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (0, 3))

if __name__ == '__main__':
    unittest.main()
