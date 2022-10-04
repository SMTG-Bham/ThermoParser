"""Tests the tp.plot.frequency module."""

import numpy as np
import unittest
import tp
from unittest.mock import Mock, patch
from tp.plot import frequency
import matplotlib.pyplot as plt

class DosTest(unittest.TestCase):
    def setUp(self):
        self.data =   {'frequency': [0, 1],
                       'A':         [0, 1],
                       'B':         [0, 1],
                       'total':     [0, 2],
                       'meta':      {}}
        self.colour = {'A':         '#ff0000',
                       'B':         '#00ff0080',
                       'total':     '#000000'}
        self.ax = Mock()

    def test_default(self):
        frequency.add_dos(self.ax, self.data, colour=self.colour, total=False,
                          main=True, scale=False, fill=False, line=True,
                          invert=False)

        self.assertEqual(self.ax.plot.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_called_once()

    def test_fill(self):
        frequency.add_dos(self.ax, self.data, colour=self.colour, total=False,
                          main=True, scale=False, fill=True, line=True,
                          invert=False)

        self.assertEqual(self.ax.fill_between.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_called_once()

    def test_not_main(self):
        frequency.add_dos(self.ax, self.data, colour=self.colour, total=False,
                          main=False, scale=False, fill=False, line=True,
                          invert=False)

        self.assertEqual(self.ax.plot.call_count, 2)
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()

    def test_total(self):
        frequency.add_dos(self.ax, self.data, colour=self.colour, total=True,
                          main=True, scale=False, fill=False, line=True,
                          invert=False)

        self.assertEqual(self.ax.plot.call_count, 3)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_called_once()

    def test_invert(self):
        frequency.add_dos(self.ax, self.data, colour=self.colour, total=False,
                          main=True, scale=False, fill=False, line=True,
                          invert=True)

        self.assertEqual(self.ax.plot.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_not_called()
        self.ax.tick_params.assert_called_once()
        self.ax.set_xlim.assert_called_once()
        self.ax.set_ylim.assert_not_called()

    @patch.object(np, 'log10')
    def test_linear_scale(self, mock_log10):
        mock_log10.side_effect = np.log10
        self.ax.get_ylim.return_value = [0, 1]
        self.ax.get_yaxis().get_scale.return_value = 'linear'
        self.ax.reset_mock()

        frequency.add_dos(self.ax, self.data, colour=self.colour, total=False,
                          main=False, scale=True, fill=False, line=True,
                          invert=False)

        self.assertEqual(self.ax.plot.call_count, 2)
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        self.ax.get_xlim.assert_called_once()
        self.ax.get_ylim.assert_called_once()
        self.ax.get_xaxis().get_scale.assert_called_once()
        self.ax.get_yaxis().get_scale.assert_called_once()
        mock_log10.assert_not_called()

    @patch.object(np, 'log10')
    def test_log_scale(self, mock_log10):
        mock_log10.return_value = [0, 1]
        self.ax.get_ylim.return_value = [1, 10]
        self.ax.get_yaxis().get_scale.return_value = 'log'
        self.ax.reset_mock()

        frequency.add_dos(self.ax, self.data, colour=self.colour, total=False,
                          main=False, scale=True, fill=False, line=True,
                          invert=False)

        self.assertEqual(self.ax.plot.call_count, 2)
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        self.ax.get_xlim.assert_called_once()
        self.ax.get_ylim.assert_called_once()
        self.ax.get_xaxis().get_scale.assert_called_once()
        self.ax.get_yaxis().get_scale.assert_called_once()
        mock_log10.assert_called_once()

class CumKappaTest(unittest.TestCase):
    def setUp(self):
        self.data = {'frequency':   [0, 1],
                     'mode_kappa':  [0, 2],
                     'temperature': [0],
                     'meta':        {'temperature': 0}}
        self.ax = Mock()

    @patch.object(tp.data.utilities, 'resolve')
    def test_default(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_cum_kappa(self.ax, self.data, direction='x', main=True,
                                scale=False, fill=False, invert=False)

        self.ax.plot.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_called_once_with(0, 2)
        mock_resolve.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_fill(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_cum_kappa(self.ax, self.data, direction='x', main=True,
                                scale=False, fill=True, invert=False)

        self.ax.fill_between.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_called_once_with(0, 2)
        mock_resolve.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_not_main(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_cum_kappa(self.ax, self.data, direction='x', main=False,
                                scale=False, fill=False, invert=False)

        self.ax.plot.assert_called_once()
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        mock_resolve.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_invert(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_cum_kappa(self.ax, self.data, direction='x', main=True,
                                scale=False, fill=False, invert=True)

        self.ax.plot.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_not_called()
        self.ax.tick_params.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 2)
        self.ax.set_ylim.assert_called_once_with(0, 1)
        mock_resolve.assert_called_once()

    @patch.object(np, 'log10')
    @patch.object(tp.data.utilities, 'resolve')
    def test_linear_scale(self, mock_resolve, mock_log10):
        mock_resolve.return_value = self.data
        mock_log10.side_effect = np.log10
        self.ax.get_ylim.return_value = [0, 1]
        self.ax.get_yaxis().get_scale.return_value = 'linear'
        self.ax.reset_mock()

        frequency.add_cum_kappa(self.ax, self.data, direction='x', main=False,
                                scale=True, fill=False, invert=False)


        self.ax.plot.assert_called_once()
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        self.ax.get_xlim.assert_called_once()
        self.ax.get_ylim.assert_called_once()
        self.ax.get_xaxis().get_scale.assert_called_once()
        self.ax.get_yaxis().get_scale.assert_called_once()
        mock_log10.assert_not_called()

    @patch.object(np, 'log10')
    @patch.object(tp.data.utilities, 'resolve')
    def test_log_scale(self, mock_resolve, mock_log10):
        mock_resolve.return_value = self.data
        mock_log10.return_value = [0, 1]
        self.ax.get_ylim.return_value = [1, 10]
        self.ax.get_yaxis().get_scale.return_value = 'log'
        self.ax.reset_mock()

        frequency.add_cum_kappa(self.ax, self.data, direction='x', main=False,
                                scale=True, fill=False, invert=False)

        self.ax.plot.assert_called_once()
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        self.ax.get_xlim.assert_called_once()
        self.ax.get_ylim.assert_called_once()
        self.ax.get_xaxis().get_scale.assert_called_once()
        self.ax.get_yaxis().get_scale.assert_called_once()
        mock_log10.assert_called_once()

class WaterfallTest(unittest.TestCase):
    def setUp(self):
        self.data = {'frequency':   [range(1000)],
                     'mode_kappa':  [range(1000)],
                     'temperature': [0],
                     'meta':        {}}
        self.ax = Mock()

    @patch.object(tp.data.utilities, 'resolve')
    def test_default(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_waterfall(self.ax, self.data, 'mode_kappa',
                                direction='x', main=True, invert=False)

        self.assertEqual(self.ax.scatter.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 999)
        self.ax.set_ylim.assert_called_once_with(11, 999)
        self.assertEqual(mock_resolve.call_count, 2)

    @patch.object(tp.data.utilities, 'resolve')
    def test_not_main(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_waterfall(self.ax, self.data, 'mode_kappa',
                                direction='x', main=False, invert=False)

        self.ax.scatter.assert_called_once()
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        mock_resolve.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_invert(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_waterfall(self.ax, self.data, 'mode_kappa',
                                direction='x', main=True, invert=True)

        self.assertEqual(self.ax.scatter.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_not_called()
        self.ax.tick_params.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(11, 999)
        self.ax.set_ylim.assert_called_once_with(0, 999)
        self.assertEqual(mock_resolve.call_count, 2)

class DensityTest(unittest.TestCase):
    def setUp(self):
        self.data = {'frequency':   [range(1000)],
                     'lifetime':    [list(range(1000))],
                     'temperature': [0],
                     'meta':        {}}
        self.data['lifetime'][0][500] = self.data['lifetime'][0][501]
        self.ax = Mock() 
 
    @patch.object(tp.data.utilities, 'resolve')
    def test_default(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_density(self.ax, self.data, 'lifetime',
                                 main=True, invert=False)

        self.assertEqual(self.ax.scatter.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 999)
        self.ax.set_ylim.assert_called_once_with(11, 999)
        self.assertEqual(mock_resolve.call_count, 2)

    @patch.object(tp.data.utilities, 'resolve')
    def test_not_main(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_density(self.ax, self.data, 'lifetime',
                                 main=False, invert=False)

        self.ax.scatter.assert_called_once()
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        mock_resolve.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_invert(self, mock_resolve):
        mock_resolve.return_value = self.data

        frequency.add_density(self.ax, self.data, 'lifetime',
                                 main=True, invert=True)

        self.assertEqual(self.ax.scatter.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_not_called()
        self.ax.tick_params.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(11, 999)
        self.ax.set_ylim.assert_called_once_with(0, 999)
        self.assertEqual(mock_resolve.call_count, 2)

class ProjectedWaterfallTest(unittest.TestCase):
    def setUp(self):
        self.data = {'frequency':   [range(1000)],
                     'mode_kappa':  [[[[[range(1000)]]]]],
                     'temperature': [0],
                     'meta':        {}}
        self.ax = Mock()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve):
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        frequency.add_projected_waterfall(self.ax, self.data, 'mode_kappa',
                                          'mode_kappa', direction='x',
                                           main=True, invert=False)

        self.assertEqual(self.ax.scatter.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 999)
        self.ax.set_ylim.assert_called_once_with(11, 999)
        mock_colourbar.assert_called_once()
        self.assertEqual(mock_resolve.call_count, 2)

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_not_main(self, mock_colourbar, mock_resolve):
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        frequency.add_projected_waterfall(self.ax, self.data, 'mode_kappa',
                                          'mode_kappa', direction='x',
                                          main=False, invert=False)

        self.ax.scatter.assert_called_once()
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        mock_colourbar.assert_called_once()
        mock_resolve.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_invert(self, mock_colourbar, mock_resolve):
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        frequency.add_projected_waterfall(self.ax, self.data, 'mode_kappa',
                                          'mode_kappa', direction='x',
                                          main=True, invert=True)

        self.assertEqual(self.ax.scatter.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_not_called()
        self.ax.tick_params.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(11, 999)
        self.ax.set_ylim.assert_called_once_with(0, 999)
        mock_colourbar.assert_called_once()
        self.assertEqual(mock_resolve.call_count, 2)

if __name__ == '__main__':
    unittest.main()