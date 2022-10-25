"""Tests the tp.plot.mfp module."""

import unittest
import tp
from unittest.mock import MagicMock, patch
from tp.plot import mfp

class CumKappaTest(unittest.TestCase):
    def setUp(self):
        self.data = {'mean_free_path': [0, 1, 2],
                     'mode_kappa':     [0, 2, 100],
                     'temperature':    [0],
                     'meta':           {'temperature': 0}}
        self.ax = MagicMock()

    @patch.object(tp.data.utilities, 'resolve')
    def test_default(self, mock_resolve):
        mock_resolve.return_value = self.data

        mfp.add_cum_kappa(self.ax, self.data, direction='x', main=True,
                          scale=False, fill=False)

        self.ax.plot.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once()
        self.ax.set_ylim.assert_called_once()
        mock_resolve.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_fill(self, mock_resolve):
        mock_resolve.return_value = self.data

        mfp.add_cum_kappa(self.ax, self.data, direction='x', main=True,
                          scale=False, fill=True)

        self.ax.fill_between.assert_called_once()
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once()
        self.ax.set_ylim.assert_called_once()
        mock_resolve.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_not_main(self, mock_resolve):
        mock_resolve.return_value = self.data

        mfp.add_cum_kappa(self.ax, self.data, direction='x', main=False,
                          scale=False, fill=False)

        self.ax.plot.assert_called_once()
        self.ax.set_xlabel.assert_not_called()
        self.ax.set_ylabel.assert_not_called()
        self.ax.set_xlim.assert_not_called()
        self.ax.set_ylim.assert_not_called()
        mock_resolve.assert_called_once()

    @patch.object(tp.data.utilities, 'resolve')
    def test_markers(self, mock_resolve):
        mock_resolve.return_value = self.data

        mfp.add_cum_kappa(self.ax, self.data, direction='x', main=True,
                          scale=False, fill=False, xmarkers=1, ymarkers=1)

        self.assertEqual(self.ax.plot.call_count, 3)
        self.ax.set_xlabel.assert_called_once()
        self.ax.set_ylabel.assert_called_once()
        self.ax.set_xlim.assert_called_once()
        self.ax.set_ylim.assert_called_once()
        mock_resolve.assert_called_once()

if __name__ == '__main__':
    unittest.main()