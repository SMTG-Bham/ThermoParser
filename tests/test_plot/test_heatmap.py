"""Tests the tp.plot.heatmap module

Failures in tp.plot.axes, tp.data.aniso, tp.data.resolve and 
tp.calculate can have knock-on effects here.
"""

import unittest
from tp.plot import heatmap
from tp.plot.axes import one_colourbar
import matplotlib as mpl
import matplotlib.pyplot as plt

class LinearHeatmapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x = [2, 4, 6]
        cls.y = [4, 8, 12]
        cls.c = [[ 2,  4,  6],
                 [ 6,  8, 10],
                 [10, 12, 14]]
        cls.fig, cls.ax = one_colourbar()
        cls.ax, cls.cbar = heatmap.add_heatmap(cls.ax, cls.x, cls.y, cls.c,
                                               xinterp=5, yinterp=5)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (2, 7))

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), (4, 14))

class LogHeatmapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x = [1e2, 1e4, 1e6]
        cls.y = [1e4, 1e8, 1e12]
        cls.c = [[1e2,  1e4,  1e6],
                  [1e6,  1e8,  1e10],
                  [1e10, 1e12, 1e14]]
        cls.fig, cls.ax = one_colourbar()
        cls.ax, cls.cbar = heatmap.add_heatmap(cls.ax, cls.x, cls.y, cls.c,
                                               xscale='log', yscale='log',
                                               cscale='log',
                                               xinterp=5, yinterp=5)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (1e2, 1e7))

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), (1e4, 1e14))

class FixedScaleHeatmapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.x = [1, 2, 3, 4]
        cls.y = [2, 4, 6, 8]
        cls.c = [[1, 2, 3],
                  [2, 4, 6],
                  [3, 6, 9]]
        cls.fig, cls.ax = one_colourbar()
        cls.ax, cls.cbar = heatmap.add_heatmap(cls.ax, cls.x, cls.y, cls.c)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xlim(self):
        self.assertEqual(self.ax.get_xlim(), (1, 4))

    def test_ylim(self):
        self.assertEqual(self.ax.get_ylim(), (2, 8))

class ZTMapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = {'temperature': [1, 2],
                 'doping':      [10, 100],
                 'zt':          [[1, 2],
                                 [4, 8]]}
        cls.fig, cls.ax = one_colourbar()
        cls.ax, cls.cbar = heatmap.add_ztmap(cls.ax, cls.d)

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xscale(self):
        self.assertEqual(self.ax.get_xscale(), 'linear')

    def test_yscale(self):
        self.assertEqual(self.ax.get_yscale(), 'log')

class CalcultedZTMapTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        c = [[[[2.]], [[4.]]],
             [[[8.]], [[16.]]]]
        s = [[[[3e6]], [[6e6]]],
             [[[9e6]], [[12e6]]]]
        etc = [[[[3.]], [[3.]]],
               [[[3.]], [[3.]]]]
        ltc = [[6.], [6.]]
        t = [2., 4.]
        d = [10., 100.]
        cls.d = {'conductivity':                    c,
                 'seebeck':                         s,
                 'electronic_thermal_conductivity': etc,
                 'lattice_thermal_conductivity':    ltc,
                 'temperature':                     t,
                 'doping':                          d,
                 'meta':                            {'units': {}}}
        cls.fig, cls.ax = one_colourbar()
        cls.ax, cls.cbar = heatmap.add_ztmap(cls.ax, cls.d, direction='x')

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xscale(self):
        self.assertEqual(self.ax.get_xscale(), 'linear')

    def test_yscale(self):
        self.assertEqual(self.ax.get_yscale(), 'log')

class TargetKLTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        c = [[[[2.]], [[4.]]],
             [[[8.]], [[16.]]]]
        s = [[[[3e6]], [[6e6]]],
             [[[9e6]], [[12e6]]]]
        zt = [[[[4]],  [[8]]],
              [[[12]], [[16]]]]
        etc = [[[[3.]], [[3.]]],
               [[[1.]], [[6.]]]]
        t = [2., 4.]
        d = [10., 100.]
        cls.d = {'conductivity':                    c,
                 'seebeck':                         s,
                 'electronic_thermal_conductivity': etc,
                 'zt':                              zt,
                 'temperature':                     t,
                 'doping':                          d,
                 'meta':                            {'units': {}}}
        cls.fig, cls.ax = one_colourbar()
        cls.ax, cls.cbar = heatmap.add_kappa_target(cls.ax, cls.d, direction='x')

    @classmethod
    def tearDownClass(cls):
        plt.close()

    def test_xscale(self):
        self.assertEqual(self.ax.get_xscale(), 'linear')

    def test_yscale(self):
        self.assertEqual(self.ax.get_yscale(), 'log')

class ColourTest(unittest.TestCase):
    def setUp(self):
        self.x = [2, 4, 6]
        self.y = [4, 8, 12]
        self.c = [[ 2,  4,  6],
                  [ 6,  8, 10],
                  [10, 12, 14]]

    def tearDown(self):
        plt.close()

    def test_colourmap(self):
        cmap = mpl.cm.get_cmap('viridis')
        self.fig, self.ax = one_colourbar()
        self.ax, self.cbar = heatmap.add_heatmap(self.ax, self.x, self.y,
                                                 self.c, xinterp=5, yinterp=5,
                                                 colour=cmap)

    def test_colour(self):
        colour = '#ff0000'
        self.fig, self.ax = one_colourbar()
        self.ax, self.cbar = heatmap.add_heatmap(self.ax, self.x, self.y,
                                                 self.c, xinterp=5, yinterp=5,
                                                 colour=colour)

    def test_colour(self):
        colours = ['#ff0000', '#00ff00']
        self.fig, self.ax = one_colourbar()
        self.ax, self.cbar = heatmap.add_heatmap(self.ax, self.x, self.y,
                                                 self.c, xinterp=5, yinterp=5,
                                                 colour=colours)

if __name__ == '__main__':
    unittest.main()
