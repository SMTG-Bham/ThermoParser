"""Tests the tp.settings module."""

import unittest
from tp import settings

class SettingsTest(unittest.TestCase):
    def test_isdict(self):
        for name, func in settings.__dict__.items():
            if name in settings.__dir__():
                self.assertIsInstance(func(), dict, '{} failed.'.format(name))
    
    def test_isint(self):
        self.assertIsInstance(settings.get_workers(), int, 'get_workers failed.')

    def test_islist(self):
        for name, func in {'style': settings.style,
                           'large_style': settings.large_style}.items():
            self.assertIsInstance(func(), list, '{} failed.'.format(name))

if __name__ == '__main__':
    unittest.main()