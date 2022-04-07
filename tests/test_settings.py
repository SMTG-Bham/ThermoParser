"""Tests the tp.settings module."""

import unittest
from tp import settings

class SettingsTest(unittest.TestCase):
    def test_isdict(self):
        for name, func in settings.__dict__.items():
            if name in settings.__dir__():
                self.assertIsInstance(func(), dict, '{} failed.'.format(name))

if __name__ == '__main__':
    unittest.main()