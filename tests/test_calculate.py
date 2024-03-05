"""Tests the tp.calculate module."""

import unittest
import numpy as np
from tp import calculate

class CumulateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Input
        cls.x = [5., 2., 7., 1., 4., 3., 6.]
        cls.y = [5., 5., 5., 5., 5., 5., 5.]
        # Correct output
        cls.x2 = [1., 2., 3., 4., 5., 6., 7.]
        cls.y2 = [5., 10., 15., 20., 25., 30., 35.]
        # Output
        cls.x, cls.y = calculate.cumulate(cls.x, cls.y)

    def test_sorting(self):
        self.assertTrue((self.x == self.x2).all())

    def test_cumulating(self):
        self.assertTrue((self.y == self.y2).all())

class GaussianTest(unittest.TestCase):
    def test_gaussian(self):
        x = [0]
        y = 1.0000 # to 4 dp
        y2 = calculate.gaussian(x, x0=0, sigma=1)
        self.assertEqual(round(y2[0], 4), y)

class LorentzianTest(unittest.TestCase):
    def test_lorentzian(self):
        x = [0]
        y = 0.6366 # to 4 dp
        y2 = calculate.lorentzian(x, x0=0, fwhm=1)
        self.assertEqual(round(y2[0], 4), y)

class LifetimeTest(unittest.TestCase):
    def test_lifetime(self):
        g = [[[0.25e-12 * np.pi**-1]]]
        t = 1.0000 # to 4 dp
        t2 = calculate.lifetime(g, use_tprc=False)
        self.assertEqual(round(t2[0][0][0], 4), t)

class MFPTest(unittest.TestCase):
    def test_mfp(self):
        g = [[[0.25e-12 * np.pi**-1]]]
        v = [[[1]]]
        l = 1.0000 # to 4 dp
        l2 = calculate.mfp(g, v, use_tprc=False)
        self.assertEqual(round(l2[0][0][0][0], 4), l)

class FDOccupationTest(unittest.TestCase):
    def test_fd_occupation(self):
        e = [0]
        ef = [0]
        t = [300]
        o = 0.5000
        o2 = calculate.fd_occupation(e, t, ef, use_tprc=False)
        self.assertEqual(round(o2[0][0][0][0], 4), o)

class BEOccupationTest(unittest.TestCase):
    def test_be_occupation(self):
        f = [20]
        t = [300]
        o = 1.5061
        o2 = calculate.be_occupation(f, t, use_tprc=False)
        self.assertEqual(round(float(o2[0]), 4), o)

class DFDDETest(unittest.TestCase):
    def test_dfdde(self):
        e = [0]
        ef = [0]
        t = [300]
        o = -9.6704
        df = calculate.dfdde(e, ef, t, use_tprc=False)
        self.assertEqual(round(df[0][0][0][0], 4), o)

class PowerFactorZTTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c =   [2.]
        cls.s =   [3e6]
        cls.ke = [3.]
        cls.kl = [6.]
        cls.k = [9.]
        cls.t =   [2.]
        cls.pf = [18.]
        cls.zt = [4.]
        cls.d = {'conductivity':                    cls.c,
                 'seebeck':                         cls.s,
                 'electronic_thermal_conductivity': cls.ke,
                 'lattice_thermal_conductivity':    cls.kl,
                 'temperature':                     cls.t,
                 'meta':                            {'units':      {},
                                                     'dimensions': {}}}

    def test_thermal_conductivity(self):
        k = calculate.thermal_conductivity(self.ke, self.kl, use_tprc=False)
        self.assertTrue((k == self.k).all())

    def test_power_factor(self):
        pf2 = calculate.power_factor(conductivity=self.c, seebeck=self.s,
                                     use_tprc=False)
        self.assertTrue((pf2 == self.pf).all())

    def test_zt(self):
        zt2 = calculate.zt(conductivity=self.c, seebeck=self.s,
                           electronic_thermal_conductivity=self.ke,
                           lattice_thermal_conductivity=self.kl,
                           temperature=self.t, use_tprc=False)
        self.assertTrue((zt2 == self.zt).all())

    def test_kl(self):
        kl2 = calculate.kl(conductivity=self.c, seebeck=self.s,
                           electronic_thermal_conductivity=self.ke,
                           zt=self.zt, temperature=self.t, use_tprc=False)
        self.assertTrue((kl2 == self.kl).all())

    def test_power_factor_fromdict(self):
        d = calculate.power_factor_fromdict(self.d, use_tprc=False)
        self.assertTrue((d['power_factor'] == self.pf).all())

    def test_zt_fromdict(self):
        d = calculate.zt_fromdict(self.d, use_tprc=False)
        self.assertTrue((d['zt'] == self.zt).all())

    def test_kl_fromdict(self):
        d = {'conductivity':                    self.c,
             'seebeck':                         self.s,
             'electronic_thermal_conductivity': self.ke,
             'zt':                              self.zt,
             'temperature':                     self.t,
             'meta':                            {'units':      {},
                                                 'dimensions': {'seebeck': ['t', 'd']}}}
        d = calculate.kl_fromdict(d, use_tprc=False)
        self.assertEqual(d['zt'], self.zt)
        self.assertEqual(d['meta']['dimensions']['lattice_thermal_conductivity'],
                         d['meta']['dimensions']['seebeck'])

if __name__ == '__main__':
    unittest.main()
