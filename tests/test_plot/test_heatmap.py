"""Tests the tp.plot.heatmap module."""

import unittest
import tp
import numpy as np
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
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('linear')
        self.ax.set_xlim.assert_called_once_with(1, 4)
        self.ax.set_ylim.assert_called_once_with(2, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')

    @patch.object(plt, 'colorbar')
    def test_log(self, mock_colourbar):
        self.x = np.power(10., self.x)
        self.y = np.power(10., self.y)
        self.c = np.power(10., self.c)
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
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('linear')
        self.ax.set_xlim.assert_called_once_with(1, 3.5)
        self.ax.set_ylim.assert_called_once_with(2, 7)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')

    @patch.object(plt, 'colorbar')
    def test_fixed_x_y(self, mock_colourbar):
        self.x.append(4)
        self.y.append(8)
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear')

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('linear')
        self.ax.set_xlim.assert_called_once_with(1, 4)
        self.ax.set_ylim.assert_called_once_with(2, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')

    @patch.object(plt, 'colorbar')
    def test_min(self, mock_colourbar):
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear', xmin=2,
                                   ymin=4, cmin=4)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('linear')
        self.ax.set_xlim.assert_called_once_with(2, 4)
        self.ax.set_ylim.assert_called_once_with(4, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')

    @patch.object(plt, 'colorbar')
    def test_max(self, mock_colourbar):
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear', xmax=2,
                                   ymax=4, cmax=4)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('linear')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(2, 6)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')

    @patch.object(plt, 'colorbar')
    def test_colourmap(self, mock_colourbar):
        try:
            cmap = mpl.cm.get_cmap('viridis')
        except AttributeError:
            cmap = mpl.colormaps['viridis']
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear',
                                   colour=cmap)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('linear')
        self.ax.set_xlim.assert_called_once_with(1, 4)
        self.ax.set_ylim.assert_called_once_with(2, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')

    @patch.object(plt, 'colorbar')
    def test_colour(self, mock_colourbar):
        colour = '#ff0000'
        cbar = heatmap.add_heatmap(self.ax, self.x, self.y, self.c,
                                   xinterp=None, yinterp=None, xscale='linear',
                                   yscale='linear', cscale='linear',
                                   colour=colour)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('linear')
        self.ax.set_xlim.assert_called_once_with(1, 4)
        self.ax.set_ylim.assert_called_once_with(2, 8)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')

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
                     'meta':                            {'units':        {'lattice_thermal_conductivity': ['']},
                                                         'dimensions':   {'lattice_thermal_conductivity': ['']},
                                                         'kappa_source': ''}}
        self.data2 = dict(self.data)
        self.data2['zt'] = [[1, 2],
                            [4, 8]]
        self.ax = Mock()

    @patch.object(tp.calculate, 'zt_fromdict')
    @patch.object(tp.calculate, 'interpolate')
    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve, mock_interpolate, mock_zt):
        mock_zt.return_value = self.data2
        mock_interpolate.return_value = (self.data, self.data)
        mock_resolve.return_value = self.data
        cbar = heatmap.add_ztmap(self.ax, self.data, kdata=self.data,
                                 xinterp=None, yinterp=None)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')
        mock_resolve.assert_called_once()
        mock_interpolate.assert_called_once()
        mock_zt.assert_called_once()

    @patch.object(tp.calculate, 'zt_fromdict')
    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_precalculated(self, mock_colourbar, mock_resolve, mock_zt):
        mock_resolve.return_value = self.data2
        cbar = heatmap.add_ztmap(self.ax, self.data2, xinterp=None,
                                 yinterp=None)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')
        mock_resolve.assert_called_once()
        mock_zt.assert_not_called()

class PFMapTest(unittest.TestCase):
    def setUp(self):
        c =   [[1, 1],
               [1, 2]]
        s =   [[1, 1],
               [1, 1]]
        t = [1, 2]
        d = [10, 100]
        self.data = {'conductivity':                    c,
                     'seebeck':                         s,
                     'temperature':                     t,
                     'doping':                          d,
                     'meta':                            {'units':      {},
                                                         'dimensions': {}}}
        self.data2 = dict(self.data)
        self.data2['power_factor'] = [[1, 1],
                                      [1, 2]]
        self.ax = Mock()

    @patch.object(tp.calculate, 'power_factor_fromdict')
    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve, mock_pf):
        mock_pf.return_value = self.data2
        mock_resolve.return_value = self.data
        cbar = heatmap.add_pfmap(self.ax, self.data, xinterp=None, yinterp=None)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')
        mock_resolve.assert_called_once()
        mock_pf.assert_called_once()

    @patch.object(tp.calculate, 'power_factor_fromdict')
    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_precalculated(self, mock_colourbar, mock_resolve, mock_pf):
        mock_resolve.return_value = self.data2
        cbar = heatmap.add_pfmap(self.ax, self.data2, xinterp=None,
                                 yinterp=None)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')
        mock_resolve.assert_called_once()
        mock_pf.assert_not_called()

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

    @patch.object(tp.calculate, 'kl_fromdict')
    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve, mock_kl):
        mock_resolve.return_value = self.data2
        mock_kl.return_value = self.data2
        cbar = heatmap.add_kappa_target(self.ax, self.data, xinterp=None,
                                        yinterp=None)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')
        mock_resolve.assert_called_once()
        mock_kl.assert_called_once()

class ZTDiffTest(unittest.TestCase):
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
                     'meta':                            {'units':        {'lattice_thermal_conductivity': ['']},
                                                         'dimensions':   {'lattice_thermal_conductivity': ['']},
                                                         'kappa_source': ''}}
        self.data2 = dict(self.data)
        self.data2['zt'] = [[1, 2],
                            [4, 8]]
        self.ax = Mock()

    @patch.object(tp.calculate, 'zt_fromdict')
    @patch.object(tp.calculate, 'interpolate')
    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve, mock_interpolate,
                     mock_zt):
        mock_zt.return_value = self.data2
        mock_interpolate.return_value = (self.data2, self.data2)
        mock_resolve.return_value = self.data
        cbar, _, l = heatmap.add_ztdiff(self.ax, self.data, self.data,
                                        kdata1=self.data, kdata2=self.data,
                                        xinterp=None, yinterp=None,
                                        label1='1', label2='2')

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')
        self.assertEqual(mock_resolve.call_count, 2)
        self.assertEqual(mock_interpolate.call_count, 4)
        self.assertEqual(mock_zt.call_count, 2)
        self.assertEqual(l, ['1', '2'])

    @patch.object(tp.calculate, 'zt_fromdict')
    @patch.object(tp.calculate, 'interpolate')
    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_precalculated(self, mock_colourbar, mock_resolve,
                           mock_interpolate, mock_zt):
        mock_interpolate.return_value = (self.data2, self.data2)
        mock_resolve.return_value = self.data2
        cbar, _, l = heatmap.add_ztdiff(self.ax, self.data2, self.data2,
                                        kdata1=self.data, kdata2=self.data,
                                        xinterp=None, yinterp=None,
                                        label1='1', label2='2')

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')
        self.assertEqual(mock_resolve.call_count, 2)
        self.assertEqual(mock_interpolate.call_count, 2)
        mock_zt.assert_not_called()
        self.assertEqual(l, ['1', '2'])

class PFDiffTest(unittest.TestCase):
    def setUp(self):
        c =   [[1, 1],
               [1, 2]]
        s =   [[1, 1],
               [1, 1]]
        t = [1, 2]
        d = [10, 100]
        self.data = {'conductivity':                    c,
                     'seebeck':                         s,
                     'temperature':                     t,
                     'doping':                          d,
                     'meta':                            {'units':        {},
                                                         'dimensions':   {}}}
        self.data2 = dict(self.data)
        self.data2['power_factor'] = [[1, 2],
                                      [4, 8]]
        self.ax = Mock()

    @patch.object(tp.calculate, 'power_factor_fromdict')
    @patch.object(tp.calculate, 'interpolate')
    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve, mock_interpolate,
                     mock_pf):
        mock_pf.return_value = self.data2
        mock_interpolate.return_value = (self.data2, self.data2)
        mock_resolve.return_value = self.data
        cbar, _, l = heatmap.add_pfdiff(self.ax, self.data, self.data,
                                        xinterp=None, yinterp=None,
                                        label1='1', label2='2')

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')
        self.assertEqual(mock_resolve.call_count, 2)
        self.assertEqual(mock_interpolate.call_count, 2)
        self.assertEqual(mock_pf.call_count, 2)
        self.assertEqual(l, ['1', '2'])

    @patch.object(tp.calculate, 'power_factor_fromdict')
    @patch.object(tp.calculate, 'interpolate')
    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_precalculated(self, mock_colourbar, mock_resolve,
                           mock_interpolate, mock_pf):
        mock_interpolate.return_value = (self.data2, self.data2)
        mock_resolve.return_value = self.data2
        cbar, _, l = heatmap.add_pfdiff(self.ax, self.data2, self.data2,
                                        kdata1=self.data, kdata2=self.data,
                                        xinterp=None, yinterp=None,
                                        label1='1', label2='2')

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xscale.assert_called_once_with('linear')
        self.ax.set_yscale.assert_called_once_with('log')
        self.ax.set_xlim.assert_called_once_with(1, 3)
        self.ax.set_ylim.assert_called_once_with(10, 1000)
        mock_colourbar.assert_called_once()
        cbar.ax.set_yscale.assert_called_once_with('linear')
        self.assertEqual(mock_resolve.call_count, 2)
        self.assertEqual(mock_interpolate.call_count, 2)
        mock_pf.assert_not_called()
        self.assertEqual(l, ['1', '2'])

if __name__ == '__main__':
    unittest.main()
