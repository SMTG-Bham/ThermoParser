"""Tests the tp.plot.frequency.add_cum_kappa function

Failures in tp.plot.axes and tp.calculate can have knock-on effects here.
"""

import unittest
from tp.plot import frequency
from tp.plot.axes import one_small_legend
import matplotlib.pyplot as plt

class BaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'frequency':   [0, 1],
                 'mode_kappa':  [[[[[0, 2]]]]],
                 'temperature': [11],
                 'meta':        {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ax = frequency.add_cum_kappa(cls.ax, cls.d, direction='x',
                                         main=True, scale=False)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (0, 1))

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), (0, 2))

class InvertTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'frequency':  [0, 1],
                 'mode_kappa':  [[[[[0, 2]]]]],
                 'temperature': [11],
                 'meta':        {}}
        cls.fig, cls.ax = one_small_legend()
        cls.fig, cls.ax = one_small_legend()
        cls.ax = frequency.add_cum_kappa(cls.ax, cls.d, direction='x',
                                         main=True, scale=False, invert=True)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (0, 2))

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), (0, 1))

class LinearScaleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'frequency':  [0, 1],
                 'mode_kappa':  [[[[[0, 2]]]]],
                 'temperature': [11],
                 'meta':        {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ylim = (1, 2)
        cls.ax.set_ylim(cls.ylim)
        cls.ax = frequency.add_cum_kappa(cls.ax, cls.d, direction='x',
                                         main=False, scale=True)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), self.ylim)

class LogScaleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'frequency':  [0, 1],
                 'mode_kappa':  [[[[[0, 2]]]]],
                 'temperature': [11],
                 'meta':        {}}
        cls.fig, cls.ax = one_small_legend()
        cls.ylim = cls.ax.get_ylim()
        cls.ax.set_yscale('log')
        cls.ylim = (10, 100)
        cls.ax.set_ylim(cls.ylim)
        cls.ax = frequency.add_cum_kappa(cls.ax, cls.d, direction='x',
                                         main=False, scale=True)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), self.ylim)

if __name__ == '__main__':
    unittest.main()
