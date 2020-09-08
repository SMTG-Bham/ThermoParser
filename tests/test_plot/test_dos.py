"""Tests the tp.plot.frequency.add_dos function

Failures in tp.plot.axes can have knock-on effects here.
"""

import unittest
from tp.plot import frequency
from tp.plot.axes import one_small_legend
import matplotlib.pyplot as plt

class NoTotalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'frequency': [0, 1],
                 'A':         [0, 3],
                 'B':         [4, 0],
                 'total':     [4, 3]}
        cls.c = {'A':         '#ff0000',
                 'B':         '#00ff0080',
                 'total':     '#000000'}
        cls.fig, cls.ax = one_small_legend()
        cls.ax = frequency.add_dos(cls.ax, cls.d, cls.c,
                                   main=True, scale=False)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), (0, 4))

    def test_legend(self):
        labels = self.ax.get_legend_handles_labels()[1]
        self.assertEqual(labels, ['A', 'B'])

class TotalTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'frequency': [0, 1],
                 'A':         [0, 3],
                 'B':         [4, 0],
                 'total':     [4, 3]}
        cls.c = {'A':         '#ff0000',
                 'B':         '#00ff0080',
                 'total':     '#000000'}
        cls.fig, cls.ax = one_small_legend()
        cls.ax = frequency.add_dos(cls.ax, cls.d, cls.c,
                                   main=True, scale=False, total=True)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), (0, 4))

    def test_legend(self):
        labels = self.ax.get_legend_handles_labels()[1]
        self.assertEqual(labels, ['Total', 'A', 'B'])

class InvertTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'frequency': [0, 1],
                 'A':         [0, 3],
                 'B':         [4, 0],
                 'total':     [4, 3]}
        cls.c = {'A':         '#ff0000',
                 'B':         '#00ff0080',
                 'total':     '#000000'}
        cls.fig, cls.ax = one_small_legend()
        cls.ax = frequency.add_dos(cls.ax, cls.d, cls.c,
                                   main=True, scale=False, invert=True)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (0, 4))

    def test_legend(self):
        labels = self.ax.get_legend_handles_labels()[1]
        self.assertEqual(labels, ['A', 'B'])

class LinearScaleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'frequency': [0, 1],
                 'A':         [0, 3],
                 'B':         [4, 0],
                 'total':     [4, 3]}
        cls.c = {'A':         '#ff0000',
                 'B':         '#00ff0080',
                 'total':     '#000000'}
        cls.fig, cls.ax = one_small_legend()
        cls.ylim = [1, 2]
        cls.ax.set_ylim(cls.ylim)
        cls.ax = frequency.add_dos(cls.ax, cls.d, cls.c,
                                   main=False, scale=True)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_ymax(self):
        self.assertEqual(self.d['B'][1], self.ylim[0])

    def test_ymax(self):
        self.assertEqual(self.d['B'][0], self.ylim[1])

    def test_legend(self):
        labels = self.ax.get_legend_handles_labels()[1]
        self.assertEqual(labels, ['A', 'B'])

class LogScaleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'frequency': [0, 1],
                 'A':         [0, 3],
                 'B':         [4, 0],
                 'total':     [4, 3]}
        cls.c = {'A':         '#ff0000',
                 'B':         '#00ff0080',
                 'total':     '#000000'}
        cls.fig, cls.ax = one_small_legend()
        cls.ylim = cls.ax.get_ylim()
        cls.ax.set_yscale('log')
        cls.ylim = [10, 100]
        cls.ax.set_ylim(cls.ylim)
        cls.ax = frequency.add_dos(cls.ax, cls.d, cls.c,
                                   main=False, scale=True)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_ymin(self):
        self.assertEqual(self.d['B'][1], self.ylim[0])

    def test_ymax(self):
        self.assertEqual(self.d['B'][0], self.ylim[1])

    def test_legend(self):
        labels = self.ax.get_legend_handles_labels()[1]
        self.assertEqual(labels, ['A', 'B'])

if __name__ == '__main__':
    unittest.main()
