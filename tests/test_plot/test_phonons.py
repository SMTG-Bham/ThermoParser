"""Tests the tp.plot.phonons module."""

import unittest
import tp
from unittest.mock import MagicMock, patch
from tp.plot import phonons
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import warnings

class DispersionTest(unittest.TestCase):
    def setUp(self):
        self.data = {'x':             [0, 1],
                     'qpoint':        [[0, 1, 2],
                                       [3, 1, 2]],
                     'frequency':     [[0, 1, 2],
                                       [0, 1, 2]],
                     'tick_position': [0, 1],
                     'tick_label':    ['A', 'B'],
                     'meta':          {}}
        self.ax = MagicMock()

    def test_default(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1

        phonons.add_dispersion(self.ax, self.data, main=True)

        self.assertEqual(self.ax.plot.call_count, 3)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_called_once_with(bottom=0)

    def test_bandmin(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1

        phonons.add_dispersion(self.ax, self.data, main=True, bandmin=1)

        self.assertEqual(self.ax.plot.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_not_called()

    def test_bandmax(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1

        phonons.add_dispersion(self.ax, self.data, main=True, bandmax=1)

        self.ax.plot.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_called_once_with(bottom=0)

    def test_not_main(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1

        phonons.add_dispersion(self.ax, self.data, main=False)

        self.assertEqual(self.ax.plot.call_count, 3)
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()

    def test_imaginary(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        self.data['frequency'][0][0] = -1

        phonons.add_dispersion(self.ax, self.data, main=True)

        self.assertEqual(self.ax.plot.call_count, 3)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_not_called()

class MultiTest(unittest.TestCase):
    def setUp(self):
        self.data = [{'x':             [0, 1],
                      'qpoint':        [[0, 1, 2],
                                        [3, 1, 2]],
                      'frequency':     [[0, 1, 2],
                                        [2, 1, 1]],
                      'tick_position': [0, 1],
                      'tick_label':    ['A', 'B'],
                      'meta':          {}},
                     {'x':             [0, 2],
                      'qpoint':        [[0, 1, 2],
                                        [3, 1, 2]],
                      'frequency':     [[0, 1, 2],
                                        [2, 1, 3]],
                      'tick_position': [0, 2],
                      'tick_label':    ['A', 'B'],
                      'meta':          {}}]
        self.ax = MagicMock()

    def test_default(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1

        phonons.add_multi(self.ax, self.data, main=True)

        self.assertEqual(self.ax.plot.call_count, 6)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_called_once_with(bottom=0)

    def test_bandmin(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1

        phonons.add_multi(self.ax, self.data, main=True, bandmin=1)

        self.assertEqual(self.ax.plot.call_count, 4)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_called_once_with(bottom=0)

    def test_bandmax(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1

        phonons.add_multi(self.ax, self.data, main=True, bandmax=1)

        self.assertEqual(self.ax.plot.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_called_once_with(bottom=0)

    def test_not_main(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1

        phonons.add_multi(self.ax, self.data, main=False)

        self.assertEqual(self.ax.plot.call_count, 6)
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()

    def test_imaginary(self):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        self.data[0]['frequency'][0][0] = -1

        phonons.add_multi(self.ax, self.data, main=True)

        self.assertEqual(self.ax.plot.call_count, 6)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 1)
        self.ax.set_ylim.assert_not_called()

class MultiColourTest(unittest.TestCase):
    def setUp(self):
        line = {'x':             [0, 1],
                'qpoint':        [[0, 1, 2],
                                  [3, 1, 2]],
                'frequency':     [[0, 1, 2],
                                  [2, 1, 1]],
                'tick_position': [0, 1],
                'tick_label':    ['A', 'B'],
                'meta':          {}}
        # three lines so a two-colour list is read as a colourmap
        self.data = [line, line, line]
        self.ax = MagicMock()

    @patch.object(phonons, 'add_dispersion')
    def test_colour_forms(self, mock_dispersion):
        x = np.linspace(0, 1, 3)
        coolwarm = plt.get_cmap('coolwarm')
        cases = [('viridis',                plt.get_cmap('viridis')(x)),
                 (coolwarm,                 coolwarm(x)),
                 (tp.plot.colour.skelton(), tp.plot.colour.skelton()(x)),
                 ('skelton',                tp.plot.colour.skelton()(x)),
                 (['red', 'blue'],          tp.plot.colour.linear('red',
                                                                  'blue')(x)),
                 ({'cmax': 'red'},          tp.plot.colour.linear('red')(x)),
                 ('red',                    ['red'] * 3),
                 (['red', 'blue', 'green'], ['red', 'blue', 'green'])]

        for colour, expected in cases:
            with self.subTest(colour=colour):
                mock_dispersion.reset_mock()
                phonons.add_multi(self.ax, self.data, colour=colour,
                                  main=False)

                self.assertEqual(mock_dispersion.call_count, 3)
                for call, e in zip(mock_dispersion.call_args_list, expected):
                    c = call.kwargs['colour']
                    if isinstance(c, list): # colourmap colours are wrapped
                        c = c[0]
                    np.testing.assert_allclose(mpl.colors.to_rgba(c),
                                               mpl.colors.to_rgba(e))

class AltDispersionTest(unittest.TestCase):
    def setUp(self):
        # warnings ignored as there are tonnes, likely due to random
        # input data to pymatgen (but filtering pymatgen warnings failed).
        warnings.simplefilter('ignore')

        import os, re
        cwd = os.getcwd()
        self.poscar = re.match('^.*ThermoParser', cwd).group() + '/tests/data/POSCAR'
        self.data =  {'qpoint':        [[0, 1, 2.1],
                                        [3.2, 0.9, 2],
                                        [0.1, 0, 0.1],
                                        [1, 1, 1],
                                        [4, 1, 0]],
                      'mode_kappa':    [[1, 2, 3],
                                        [2, 1, 3],
                                        [1, 2, 1],
                                        [2, 1, 3],
                                        [0, 3, 4]],
                      'temperature':   [11],
                      'meta':          {}}
        self.pdata = {'x':             [0, 1, 2, 3],
                      'qpoint':        [[0, 1, 2],
                                        [3, 1, 2],
                                        [0, 0, 0],
                                        [1, 1, 1]],
                      'frequency':     [[0, 1, 2],
                                        [2, 1, 0],
                                        [3, 2, 1],
                                        [1, 2, 3]],
                      'tick_position': [0, 3],
                      'tick_label':    ['A', 'B'],
                      'meta':          {}}
        self.ax = MagicMock()

    @patch.object(tp.data.utilities, 'resolve')
    def test_default(self, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_resolve.return_value = self.data
        
        phonons.add_alt_dispersion(self.ax, self.data, self.pdata,
                                   'mode_kappa', main=True, smoothing=1,
                                   poscar=self.poscar)

        self.assertEqual(self.ax.plot.call_count, 3)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 3)
        self.ax.set_ylim.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_bandmin(self, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_resolve.return_value = self.data

        phonons.add_alt_dispersion(self.ax, self.data, self.pdata,
                                   'mode_kappa', main=True, smoothing=1,
                                   poscar=self.poscar, bandmin=1)

        self.assertEqual(self.ax.plot.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 3)
        self.ax.set_ylim.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_bandmax(self, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_resolve.return_value = self.data

        phonons.add_alt_dispersion(self.ax, self.data, self.pdata,
                                   'mode_kappa', main=True, smoothing=1,
                                   poscar=self.poscar, bandmax=1)

        self.ax.plot.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 3)
        self.ax.set_ylim.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_not_main(self, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_resolve.return_value = self.data

        phonons.add_alt_dispersion(self.ax, self.data, self.pdata,
                                   'mode_kappa', main=False, smoothing=1,
                                   poscar=self.poscar)

        self.assertEqual(self.ax.plot.call_count, 3)
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()

class ProjectedDispersionTest(unittest.TestCase):
    def setUp(self):
        # warnings ignored as there are tonnes, likely due to random
        # input data to pymatgen (but filtering pymatgen warnings failed).
        warnings.simplefilter('ignore')
        import os, re
        cwd = os.getcwd()
        self.poscar = re.match('^.*ThermoParser', cwd).group() + '/tests/data/POSCAR'
        self.data =  {'qpoint':        [[0, 1, 2.1],
                                        [3.2, 0.9, 2],
                                        [0.1, 0, 0.1],
                                        [1, 1, 1],
                                        [4, 1, 0]],
                      'mode_kappa':    [[1, 2, 3],
                                        [2, 1, 3],
                                        [1, 2, 1],
                                        [2, 1, 3],
                                        [0, 3, 4]],
                      'temperature':   [11],
                      'meta':          {}}
        self.pdata = {'x':             [0, 1, 2, 3, 4, 5],
                      'qpoint':        [[0, 1, 2],
                                        [3, 1, 2],
                                        [0, 0, 0],
                                        [3, 1, 2],
                                        [0, 0, 0],
                                        [1, 1, 1]],
                      'frequency':     [[0, 1, 2],
                                        [2, 1, 0],
                                        [3, 2, 1],
                                        [2, 1, 0],
                                        [3, 2, 1],
                                        [1, 2, 3]],
                      'tick_position': [0, 5],
                      'tick_label':    ['A', 'B'],
                      'meta':          {}}
        self.ax = MagicMock()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        phonons.add_projected_dispersion(self.ax, self.data, self.pdata,
                                         'mode_kappa', main=True, smoothing=1,
                                         poscar=self.poscar)

        self.assertEqual(self.ax.scatter.call_count, 3)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 5)
        self.ax.set_ylim.assert_called_once()
        mock_colourbar.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_bandmin(self, mock_colourbar, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        phonons.add_projected_dispersion(self.ax, self.data, self.pdata,
                                         'mode_kappa', main=True, smoothing=1,
                                         poscar=self.poscar, bandmin=1)

        self.assertEqual(self.ax.scatter.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 5)
        self.ax.set_ylim.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_bandmax(self, mock_colourbar, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        phonons.add_projected_dispersion(self.ax, self.data, self.pdata,
                                         'mode_kappa', main=True, smoothing=1,
                                         poscar=self.poscar, bandmax=1)

        self.ax.scatter.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 5)
        self.ax.set_ylim.assert_called_once()
        mock_colourbar.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_not_main(self, mock_colourbar, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        phonons.add_projected_dispersion(self.ax, self.data, self.pdata,
                                         'mode_kappa', main=False, smoothing=1,
                                         poscar=self.poscar)

        self.assertEqual(self.ax.scatter.call_count, 3)
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        mock_colourbar.assert_called_once()

class AltProjectedDispersionTest(unittest.TestCase):
    def setUp(self):
        # warnings ignored as there are tonnes, likely due to random
        # input data to pymatgen (but filtering pymatgen warnings failed).
        warnings.simplefilter('ignore')
        import os, re
        cwd = os.getcwd()
        self.poscar = re.match('^.*ThermoParser', cwd).group() + '/tests/data/POSCAR'
        self.data =  {'qpoint':        [[0, 1, 2.1],
                                        [3.2, 0.9, 2],
                                        [0.1, 0, 0.1],
                                        [1, 1, 1],
                                        [4, 1, 0]],
                      'mode_kappa':    [[1, 2, 3],
                                        [2, 1, 3],
                                        [1, 2, 1],
                                        [2, 1, 3],
                                        [0, 3, 4]],
                      'temperature':   [11],
                      'meta':          {}}
        self.pdata = {'x':             [0, 1, 2, 3, 4, 5],
                      'qpoint':        [[0, 1, 2],
                                        [3, 1, 2],
                                        [0, 0, 0],
                                        [3, 1, 2],
                                        [0, 0, 0],
                                        [1, 1, 1]],
                      'frequency':     [[0, 1, 2],
                                        [2, 1, 0],
                                        [3, 2, 1],
                                        [2, 1, 0],
                                        [3, 2, 1],
                                        [1, 2, 3]],
                      'tick_position': [0, 5],
                      'tick_label':    ['A', 'B'],
                      'meta':          {}}
        self.ax = MagicMock()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_default(self, mock_colourbar, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        phonons.add_alt_projected_dispersion(self.ax, self.data, self.pdata,
                                             'mode_kappa', 'mode_kappa',
                                             main=True, smoothing=1,
                                             poscar=self.poscar)

        self.assertEqual(self.ax.scatter.call_count, 3)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 5)
        self.ax.set_ylim.assert_called_once()
        mock_colourbar.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_bandmin(self, mock_colourbar, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        phonons.add_alt_projected_dispersion(self.ax, self.data, self.pdata,
                                             'mode_kappa', 'mode_kappa',
                                             main=True, smoothing=1,
                                             poscar=self.poscar, bandmin=1)

        self.assertEqual(self.ax.scatter.call_count, 2)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 5)
        self.ax.set_ylim.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_bandmax(self, mock_colourbar, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        phonons.add_alt_projected_dispersion(self.ax, self.data, self.pdata,
                                             'mode_kappa', 'mode_kappa',
                                             main=True, smoothing=1,
                                             poscar=self.poscar, bandmax=1)

        self.ax.scatter.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 5)
        self.ax.set_ylim.assert_called_once()
        mock_colourbar.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    @patch.object(plt, 'colorbar')
    def test_not_main(self, mock_colourbar, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_colourbar().ax.yaxis.get_scale.return_value = 'log'
        mock_colourbar.reset_mock()
        mock_resolve.return_value = self.data

        phonons.add_alt_projected_dispersion(self.ax, self.data, self.pdata,
                                             'mode_kappa', 'mode_kappa',
                                             main=False, smoothing=1,
                                             poscar=self.poscar)

        self.assertEqual(self.ax.scatter.call_count, 3)
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        mock_colourbar.assert_called_once()

class WidebandTest(unittest.TestCase):
    def setUp(self):
        # warnings ignored as there are tonnes, likely due to random
        # input data to pymatgen (but filtering pymatgen warnings failed).
        warnings.simplefilter('ignore')
        import os, re
        cwd = os.getcwd()
        self.poscar = re.match('^.*ThermoParser', cwd).group() + '/tests/data/POSCAR'
        self.data =  {'qpoint':        [[0, 1, 2.1],
                                        [3.2, 0.9, 2],
                                        [0.1, 0, 0.1],
                                        [1, 1, 1],
                                        [4, 1, 0]],
                      'gamma':         [[1, 2, 3],
                                        [2, 1, 3],
                                        [1, 2, 1],
                                        [2, 1, 3],
                                        [0, 3, 4]],
                      'temperature':   [11],
                      'meta':          {}}
        self.pdata = {'x':             [0, 1, 2, 3, 4, 5],
                      'qpoint':        [[0, 1, 2],
                                        [3, 1, 2],
                                        [0, 0, 0],
                                        [3, 1, 2],
                                        [0, 0, 0],
                                        [1, 1, 1]],
                      'frequency':     [[0, 1, 2],
                                        [2, 1, 0],
                                        [3, 2, 1],
                                        [2, 1, 0],
                                        [3, 2, 1],
                                        [1, 2, 3]],
                      'tick_position': [0, 5],
                      'tick_label':    ['A', 'B'],
                      'meta':          {}}
        self.ax = MagicMock()

    @patch.object(tp.data.utilities, 'resolve')
    def test_default(self, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_resolve.return_value = self.data

        phonons.add_wideband(self.ax, self.data, self.pdata, smoothing=1,
                             main=True, poscar=self.poscar)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once_with(0, 5)

    @patch.object(tp.data.utilities, 'resolve')
    def test_not_main(self, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_resolve.return_value = self.data

        phonons.add_wideband(self.ax, self.data, self.pdata, smoothing=1,
                             main=False, poscar=self.poscar)

        self.ax.pcolormesh.assert_called_once()
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()

    @patch.object(tp.data.utilities, 'resolve')
    def test_colour_forms(self, mock_resolve):
        self.ax.spines['bottom'].get_linewidth().return_value = 1
        mock_resolve.return_value = self.data
        coolwarm = plt.get_cmap('coolwarm')

        for colour in ['viridis', coolwarm, tp.plot.colour.skelton(),
                       '#ff0000', ['red', 'blue'], {'cmax': 'red'}]:
            with self.subTest(colour=colour):
                self.ax.reset_mock()
                phonons.add_wideband(self.ax, self.data, self.pdata,
                                     smoothing=1, main=False,
                                     poscar=self.poscar, colour=colour)

                self.ax.pcolormesh.assert_called_once()
                cmap = self.ax.pcolormesh.call_args.kwargs['cmap']
                self.assertIsInstance(cmap, mpl.colors.Colormap)
                if not isinstance(colour, str):
                    continue
                if colour == 'viridis':
                    self.assertEqual(cmap.name, 'viridis')
                else:
                    np.testing.assert_allclose(cmap(1.),
                                               mpl.colors.to_rgba(colour))

class GetEquivalentQpointTest(unittest.TestCase):
    @patch.object(np, 'dot')
    def test_get_equivalent_qpoint(self, mock_dot):
        qk = [[0, 1, 2.1],
              [3.2, 0.9, 2],
              [0.1, 0, 0.1],
              [1, 1, 1],
              [4, 1, 0]]
        qp = [[0, 1, 2],
              [3, 1, 2],
              [0, 0, 0],
              [1, 1, 1]]
        symops = MagicMock()
        def same(x, *args, **kwargs): return [x]
        mock_dot.side_effect = same

        i = [phonons.get_equivalent_qpoint(qk, symops, q) for q in qp]

        self.assertEqual(i, [0, 1, 2, 3])
        self.assertEqual(mock_dot.call_count, 4)

class TilePropertiesTest(unittest.TestCase):
    def test_single(self):
        p = 'a'
        ps = phonons.tile_properties(p, 0, 3)
        self.assertTrue((ps == ['a', 'a', 'a']).all())

    def test_single_array(self):
        p = ['a']
        ps = phonons.tile_properties(p, 0, 3)
        self.assertTrue((ps == ['a', 'a', 'a']).all())

    def test_full(self):
        p = ['a', 'b', 'c']
        ps = phonons.tile_properties(p, 0, 3)
        self.assertEqual(ps, p)

    def test_fill(self):
        p = ['a', 'b']
        ps = phonons.tile_properties(p, 0, 3)
        self.assertEqual(ps, ['a', 'b', 'b'])

    def test_subset(self):
        p = ['a', 'b', 'c']
        ps = phonons.tile_properties(p, 1, 3)
        self.assertEqual(ps, ['b', 'c'])

    def test_not_subset(self):
        p = ['a', 'b', 'c']
        ps = phonons.tile_properties(p, 2, 4)
        self.assertEqual(ps, ['a', 'b'])

if __name__ == '__main__':
    unittest.main()