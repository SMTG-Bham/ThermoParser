"""Tests the tp.plot.heatmap module."""

import unittest
import tp
import numpy as np
from glob import glob
from os import remove
from unittest.mock import Mock, patch
from tp.plot import heatmap
import matplotlib as mpl
import matplotlib.pyplot as plt

class HeatmapTest(unittest.TestCase):
    def setUp(self):
        self.x = [1, 2, 3]
        self.y = [2, 4, 6]
        self.c = [[1, 2, 3],
                  [3, 4, 5],
                  [5, 6, 7]]
        self.ax = Mock()

    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar):
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear')

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_not_called()
        self.ax.set_xlim.assert_called_once_with(1, 4)
        self.ax.set_ylim.assert_called_once_with(2, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()

    @patch.object(plt, 'colorbar')
    def test_log(self, mock_colourbar):
        self.x = np.power(10, self.x)
        self.y = np.power(10, self.y)
        self.c = np.power(10, self.c)
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='log',
                                   yscale='log', cscale='log')

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('log')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1e1, 1e4)
        self.ax.set_ylim.assert_called_once_with(1e2, 1e8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('log')

    @patch.object(plt, 'colorbar')
    def test_interpolate(self, mock_colourbar):
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=5, yinterp=5, kind='linear',
                                   xscale='linear', yscale='linear',
                                   cscale='linear')

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_not_called()
        self.ax.set_xlim.assert_called_once_with(1, 3.5)
        self.ax.set_ylim.assert_called_once_with(2, 7)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()

    @patch.object(plt, 'colorbar')
    def test_fixed_x_y(self, mock_colourbar):
        self.x.append(4)
        self.y.append(8)
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear')

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_not_called()
        self.ax.set_xlim.assert_called_once_with(1, 4)
        self.ax.set_ylim.assert_called_once_with(2, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()

    @patch.object(plt, 'colorbar')
    def test_min(self, mock_colourbar):
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear', xmin=2,
                                   ymin=4, cmin=4)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_not_called()
        self.ax.set_xlim.assert_called_once_with(2, 4)
        self.ax.set_ylim.assert_called_once_with(4, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()

    @patch.object(plt, 'colorbar')
    def test_max(self, mock_colourbar):
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear', xmax=2,
                                   ymax=4, cmax=4)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_not_called()
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(2, 6)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()

    @patch.object(plt, 'colorbar')
    def test_colourmap(self, mock_colourbar):
        cmap = mpl.cm.get_cmap('viridis')
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear',
                                   colour=cmap)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_not_called()
        self.ax.set_xlim.assert_called_once_with(1, 4)
        self.ax.set_ylim.assert_called_once_with(2, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()

    @patch.object(plt, 'colorbar')
    def test_colour(self, mock_colourbar):
        colour = '#ff0000'
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear',
                                   colour=colour)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_not_called()
        self.ax.set_xlim.assert_called_once_with(1, 4)
        self.ax.set_ylim.assert_called_once_with(2, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()

    @patch.object(plt, 'colorbar')
    def test_colours(self, mock_colourbar):
        colours = ['#ff0000', '#00ff00']
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear',
                                   colour=colours)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_not_called()
        self.ax.set_xlim.assert_called_once_with(1, 4)
        self.ax.set_ylim.assert_called_once_with(2, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()

class ZTMapTest(unittest.TestCase):
    def setUp(self):
        c =   [[1, 1],
               [1, 2]]
        s =   [[1, 1],
               [1, 1]]
        etc = [[1, 1],
               [1, 1]]
        ltc =  [1, 1]
        t = [1, 2]
        d = [10, 100]
        self.data = {'conductivity':                    c,
                     'seebeck':                         s,
                     'electronic_thermal_conductivity': etc,
                     'lattice_thermal_conductivity':    ltc,
                     'temperature':                     t,
                     'doping':                          d,
                     'meta':                            {'units': {}}}
        self.data2 = dict(self.data)
        self.data2['zt'] = [[1, 2],
                            [4, 8]]
        self.ax = Mock()

    def tearDown(self):
        dat = glob("zt.hdf5")
        for f in dat: remove(f)

    @patch.object(tp.calculate, 'zt_fromdict')
    @patch.object(tp.data.resolve, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve, mock_zt):
        mock_resolve.return_value = self.data
        mock_zt.return_value = self.data2
        cbar = heatmap.add_ztmap(self.ax, self.data, xinterp=None,
                                 yinterp=None)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()
        mock_resolve.assert_called_once()
        mock_zt.assert_called_once()

    @patch.object(tp.calculate, 'zt_fromdict')
    @patch.object(tp.data.resolve, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_precalculated(self, mock_colourbar, mock_resolve, mock_zt):
        mock_resolve.return_value = self.data2
        cbar = heatmap.add_ztmap(self.ax, self.data2, xinterp=None,
                                 yinterp=None)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()
        mock_resolve.assert_not_called()
        mock_zt.assert_not_called()

class TargetKLTest(unittest.TestCase):
    def setUp(self):
        c =   [[1, 1],
               [1, 2]]
        s =   [[1, 1],
               [1, 1]]
        etc = [[1, 1],
               [1, 1]]
        t = [1, 2]
        d = [10, 100]
        self.data = {'conductivity':                    c,
                     'seebeck':                         s,
                     'electronic_thermal_conductivity': etc,
                     'temperature':                     t,
                     'doping':                          d,
                     'meta':                            {'units': {}}}
        self.data2 = dict(self.data)
        self.data2['lattice_thermal_conductivity'] = [[1, 2],
                                                      [4, 8]]
        self.ax = Mock()

    def tearDown(self):
        dat = glob("target-kl.hdf5")
        for f in dat: remove(f)

    @patch.object(tp.calculate, 'kl_fromdict')
    @patch.object(tp.data.resolve, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve, mock_kl):
        mock_resolve.return_value = self.data
        mock_kl.return_value = self.data2
        cbar = heatmap.add_kappa_target(self.ax, self.data, xinterp=None,
                                        yinterp=None)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_not_called()
        self.ax.set_yscale.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_not_called()
        mock_resolve.assert_called_once()
        mock_kl.assert_called_once()

if __name__ == '__main__':
    unittest.main()
