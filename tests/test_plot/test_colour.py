"""Tests the tp.plot.colour module."""

import unittest
from tp.plot import colour
import matplotlib

class LinearTest(unittest.TestCase):
    def test_is_colourmap(self):
        cmap = colour.linear('#ff0000')
        self.assertIsInstance(cmap, matplotlib.colors.ListedColormap)

    def test_midpoint_colour(self):
        c = colour.linear('#000000')(0.5)
        cround = [round(c[0], 2), round(c[1], 2), round(c[2], 2)]
        self.assertEqual(cround, [0.5, 0.5, 0.5])

class ElbowTest(unittest.TestCase):
    def test_is_colourmap(self):
        cmap = colour.elbow('#ff0000')
        self.assertIsInstance(cmap, matplotlib.colors.ListedColormap)

class HighlightTest(unittest.TestCase):
    def test_is_colourmap(self):
        cmap = colour.highlight(matplotlib.pyplot.get_cmap('viridis'), 'red')
        self.assertIsInstance(cmap, matplotlib.colors.ListedColormap)

class SkeltonTest(unittest.TestCase):
    def test_is_colourmap(self):
        c = colour.skelton()
        self.assertIsInstance(c, matplotlib.colors.ListedColormap)

if __name__ == '__main__':
    unittest.main()