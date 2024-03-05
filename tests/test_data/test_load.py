"""Tests the tp.data.load module."""

import h5py
import json
import unittest
import numpy as np
import yaml
from os import remove
from pymatgen.io.vasp.inputs import Poscar
from unittest.mock import call, MagicMock, mock_open, patch
from tp.data import load
from tp import settings

class AmsetTest(unittest.TestCase):
    @patch.object(settings, 'amset_conversions')
    @patch.object(settings, 'conversions')
    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_conductivity(self, mock_opens, mock_json, mock_conversions,
                          mock_aconversions):
        mock_aconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'conductivity'
        data = {q:              np.zeros((2, 2, 3, 3)),
                'temperatures': [0, 1],
                'doping':       [1, 2]}
        data[q][0,1,0,0] = 1
        mock_json.return_value = data

        data2 = load.amset('mock', q)
        mock_opens.assert_called_once()
        mock_json.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)
        self.assertEqual(data2[q][1,0,0,0], 1)

    @patch.object(settings, 'amset_conversions')
    @patch.object(settings, 'conversions')
    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_seebeck(self, mock_opens, mock_json, mock_conversions,
                     mock_aconversions):
        mock_aconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'seebeck'
        data = {q:              np.zeros((2, 2, 3, 3)),
                'temperatures': [0, 1],
                'doping':       [1, 2]}
        data[q][0,1,0,0] = 1
        mock_json.return_value = data

        data2 = load.amset('mock', q)
        mock_opens.assert_called_once()
        mock_json.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)
        self.assertEqual(data2[q][1,0,0,0], 1)

    @patch.object(settings, 'amset_conversions')
    @patch.object(settings, 'conversions')
    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_electronic_thermal_conductivity(self, mock_opens, mock_json,
                                             mock_conversions,
                                             mock_aconversions):
        mock_aconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'electronic_thermal_conductivity'
        data = {q:              np.zeros((2, 2, 3, 3)),
                'temperatures': [0, 1],
                'doping':       [1, 2]}
        data[q][0,1,0,0] = 1
        mock_json.return_value = data

        data2 = load.amset('mock', q)
        mock_opens.assert_called_once()
        mock_json.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)
        self.assertEqual(data2[q][1,0,0,0], 1)

    @patch.object(settings, 'amset_conversions')
    @patch.object(settings, 'conversions')
    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_fermi_levels(self, mock_opens, mock_json, mock_conversions,
                          mock_aconversions):
        mock_aconversions.return_value = {}
        mock_conversions.return_value = {}
        aq = 'fermi_levels'
        tq = 'fermi_level'
        data = {aq:             np.zeros((2, 2)),
                'temperatures': [0, 1],
                'doping':       [1, 2]}
        data[aq][0,1] = 1
        mock_json.return_value = data

        data2 = load.amset('mock', aq)
        mock_opens.assert_called_once()
        mock_json.assert_called_once()
        for q2 in [tq, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)
        self.assertEqual(data2[tq][1,0], 1)

    @patch.object(settings, 'amset_conversions')
    @patch.object(settings, 'conversions')
    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_mobility(self, mock_opens, mock_json, mock_conversions,
                      mock_aconversions):
        mock_aconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'mobility'
        data = {q:              {'test':    np.zeros((2, 2, 3, 3)),
                                 'overall': np.zeros((2, 2, 3, 3))},
                'temperatures': [0, 1],
                'doping':       [1, 2]}
        data[q]['test'][0,1,0,0] = 1
        data[q]['overall'][0,1,0,0] = 1
        mock_json.return_value = data

        data2 = load.amset('mock', q)
        mock_opens.assert_called_once()
        mock_json.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)
        self.assertEqual(np.array(data2[q])[0,1,0,0,0], 1)
        self.assertEqual(list(data2['stype']), ['test', 'Total'])

class AmsetMeshTest(unittest.TestCase):

    def tearDown(self):
        remove('test.hdf5')

    @patch.object(settings, 'amset_conversions')
    @patch.object(settings, 'conversions')
    def test_scattering(self, mock_conversions, mock_aconversions):
        mock_aconversions.return_value = {}
        mock_conversions.return_value = {}
        with h5py.File('test.hdf5', 'w') as f:
            f['scattering_rates_up']              = np.zeros((1, 2, 2, 3, 1))
            f['scattering_rates_down']            = np.zeros((1, 2, 2, 3, 1))
            f['scattering_labels']                = ['test'.encode('ascii', 'ignore')]
            f['temperatures']                     = [0, 1]
            f['doping']                           = [1, 2]
            f['ir_kpoints']                       = [[1, 2, 3]]
            f['scattering_rates_up'][0,0,1,0,0]   = 1
            f['scattering_rates_down'][0,0,1,0,0] = 3

        data = load.amset_mesh('test.hdf5', 'scattering_rates', spin='up')
        for q2 in ['scattering_rates', 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data)
        self.assertEqual(data['scattering_rates'][0,1,0,0,0], 1)
        self.assertEqual(data['stype'], ['test', 'Total'])

    @patch.object(settings, 'amset_conversions')
    @patch.object(settings, 'conversions')
    def test_spin(self, mock_conversions, mock_aconversions):
        mock_aconversions.return_value = {}
        mock_conversions.return_value = {}
        with h5py.File('test.hdf5', 'w') as f:
            f['scattering_rates_up']              = np.zeros((1, 2, 2, 3, 1))
            f['scattering_rates_down']            = np.zeros((1, 2, 2, 3, 1))
            f['scattering_labels']                = ['test'.encode('ascii', 'ignore')]
            f['temperatures']                     = [0, 1]
            f['doping']                           = [1, 2]
            f['ir_kpoints']                       = [[1, 2, 3]]
            f['scattering_rates_up'][0,0,1,0,0]   = 1
            f['scattering_rates_down'][0,0,1,0,0] = 3
            f['energies_up']                      = np.zeros((1, 2, 2, 3, 1))
            f['energies_down']                    = np.zeros((1, 2, 2, 3, 1))

        data = load.amset_mesh('test.hdf5', 'scattering_rates', spin='avg')
        for q2 in ['scattering_rates', 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data)
        self.assertEqual(data['scattering_rates'][0,1,0,0,0], 2)
        self.assertEqual(data['stype'], ['test', 'Total'])

class BoltzTraPTest(unittest.TestCase):

    def tearDown(self):
        remove('test.hdf5')

    @patch.object(settings, 'boltztrap_conversions')
    @patch.object(settings, 'conversions')
    def test_conductivity(self, mock_conversions, mock_bconversions):
        mock_bconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'conductivity'
        with h5py.File('test.hdf5', 'w') as f:
            g = f.create_group(q)
            g['n']        = np.zeros((2, 2, 3, 3))
            f['temperature'] = [0, 1]
            f['doping']      = [1, 2]

        data = load.boltztrap('test.hdf5', q)
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'boltztrap_conversions')
    @patch.object(settings, 'conversions')
    def test_effective_mass(self, mock_conversions, mock_bconversions):
        mock_bconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'average_eff_mass'
        with h5py.File('test.hdf5', 'w') as f:
            g = f.create_group(q)
            g['n']           = np.zeros((2, 2, 3, 3))
            f['temperature'] = [0, 1]
            f['doping']      = [1, 2]

        data = load.boltztrap('test.hdf5', q)
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'boltztrap_conversions')
    @patch.object(settings, 'conversions')
    def test_fermi_level(self, mock_conversions, mock_bconversions):
        mock_bconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'fermi_level'
        with h5py.File('test.hdf5', 'w') as f:
            g = f.create_group(q)
            g['n']        = np.zeros((2, 2, 3, 3))
            f['temperature'] = [0, 1]
            f['doping']      = [1, 2]

        data = load.boltztrap('test.hdf5', q)
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'boltztrap_conversions')
    @patch.object(settings, 'conversions')
    def test_power_factor(self, mock_conversions, mock_bconversions):
        mock_bconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'power_factor'
        with h5py.File('test.hdf5', 'w') as f:
            g = f.create_group(q)
            g['n']        = np.zeros((2, 2, 3, 3))
            f['temperature'] = [0, 1]
            f['doping']      = [1, 2]

        data = load.boltztrap('test.hdf5', q)
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'boltztrap_conversions')
    @patch.object(settings, 'conversions')
    def test_seebeck(self, mock_conversions, mock_bconversions):
        mock_bconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'seebeck'
        with h5py.File('test.hdf5', 'w') as f:
            g = f.create_group(q)
            g['n']        = np.zeros((2, 2, 3, 3))
            f['temperature'] = [0, 1]
            f['doping']      = [1, 2]

        data = load.boltztrap('test.hdf5', q)
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'boltztrap_conversions')
    @patch.object(settings, 'conversions')
    def test_thermal_conductivity(self, mock_conversions, mock_bconversions):
        mock_bconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'electronic_thermal_conductivity'
        with h5py.File('test.hdf5', 'w') as f:
            g = f.create_group(q)
            g['n']        = np.zeros((2, 2, 3, 3))
            f['temperature'] = [0, 1]
            f['doping']      = [1, 2]

        data = load.boltztrap('test.hdf5', q)
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data)

class Phono3pyTest(unittest.TestCase):

    def tearDown(self):
        remove('test.hdf5')

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_frequency(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'frequency'
        with h5py.File('test.hdf5', 'w') as f:
            f[q]             = np.zeros((2, 3))
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in [q, 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_gamma(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'gamma'
        with h5py.File('test.hdf5', 'w') as f:
            f[q]             = np.zeros((2, 3, 3))
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in [q, 'temperature', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_group_velocity(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'group_velocity'
        with h5py.File('test.hdf5', 'w') as f:
            f[q]             = np.zeros((2, 3, 3))
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in [q, 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_gv_by_gv(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'gv_by_gv'
        with h5py.File('test.hdf5', 'w') as f:
            f[q]             = np.zeros((2, 3, 6))
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in [q, 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_heat_capacity(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'heat_capacity'
        with h5py.File('test.hdf5', 'w') as f:
            f[q]             = np.zeros((2, 3, 3))
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in [q, 'temperature', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_kappa(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'kappa'
        with h5py.File('test.hdf5', 'w') as f:
            f[q]             = np.zeros((2, 6))
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in ['lattice_thermal_conductivity', 'temperature', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_mode_kappa(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'mode_kappa'
        with h5py.File('test.hdf5', 'w') as f:
            f[q]             = np.ones((2, 2, 3, 6))
            f['kappa']       = np.full((2, 6), 6)
            f['weight']      = np.ones(2)
            f['mesh']        = np.array([2, 1, 1])
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in [q, 'temperature', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_old_mode_kappa(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        # in old phono3py versions, mode_kappa is multiplied by the mesh
        q = 'mode_kappa'
        with h5py.File('test.hdf5', 'w') as f:
            f[q]             = np.ones((2, 2, 3, 6))
            f['kappa']       = np.full((2, 6), 3)
            f['weight']      = np.ones(2)
            f['mesh']        = np.array([2, 1, 1])
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in [q, 'temperature', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_ph_ph_strength(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'ave_pp'
        with h5py.File('test.hdf5', 'w') as f:
            f[q]             = np.zeros((2, 3))
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in ['ph_ph_strength', 'meta']:
            self.assertIn(q2, data)

    @patch.object(settings, 'phono3py_conversions')
    @patch.object(settings, 'conversions')
    def test_wideband(self, mock_conversions, mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        q = 'wideband'
        with h5py.File('test.hdf5', 'w') as f:
            f['frequency']   = np.zeros((2, 3))
            f['gamma']       = np.zeros((2, 2, 3))
            f['temperature'] = np.array([0, 1])
            f['qpoint']      = np.array([0, 1])

        data = load.phono3py('test.hdf5', q)
        for q2 in ['frequency', 'gamma', 'temperature', 'qpoint', 'meta']:
            self.assertIn(q2, data)

class PhonopyDispersionTest(unittest.TestCase):
    @patch.object(settings, 'phonopy_conversions')
    @patch.object(settings, 'conversions')
    @patch.object(yaml, 'safe_load')
    @patch("builtins.open", new_callable=mock_open)
    def test_default(self, mock_opens, mock_yaml, mock_conversions,
                     mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        mock_data = MagicMock()
        mock_yaml.return_value = mock_data

        data2 = load.phonopy_dispersion('mock')
        mock_opens.assert_called_once()
        mock_yaml.assert_called_once()
        self.assertIn(call('distance'),
                      mock_data['phonon'][0].__getitem__.call_args_list)

    @patch.object(settings, 'phonopy_conversions')
    @patch.object(settings, 'conversions')
    @patch.object(yaml, 'safe_load')
    @patch("builtins.open", new_callable=mock_open)
    def test_scale(self, mock_opens, mock_yaml, mock_conversions,
                   mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        mock_data = MagicMock()
        mock_data2 = MagicMock()
        mock_yaml.return_value = mock_data

        data2 = load.phonopy_dispersion('mock', mock_data2)
        mock_opens.assert_called_once()
        mock_yaml.assert_called_once()
        self.assertIn(call('distance'),
                      mock_data['phonon'][0].__getitem__.call_args_list)
        self.assertIn(call('tick_position'),
                      mock_data2.__getitem__.call_args_list)

class PhonopyDoSTest(unittest.TestCase):
    @patch.object(settings, 'phonopy_conversions')
    @patch.object(settings, 'conversions')
    @patch.object(Poscar, 'from_file')
    @patch.object(np, 'loadtxt')
    def test_default(self, mock_load, mock_poscar, mock_conversions,
                     mock_pconversions):
        mock_pconversions.return_value = {}
        mock_conversions.return_value = {}
        mock_data = MagicMock()
        mock_load.return_value = mock_data
        mock_POSCAR = MagicMock()
        mock_poscar.return_value = mock_POSCAR

        data2 = load.phonopy_dos('mock')
        mock_load.assert_called_once()
        mock_poscar.assert_called_once()
