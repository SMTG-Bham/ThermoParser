"""Tests the tp.data.load module."""

import h5py
import json
import unittest
import numpy as np
import yaml
from pymatgen.io.vasp.inputs import Poscar
from unittest.mock import call, MagicMock, mock_open, patch
from tp.data import load

class AmsetTest(unittest.TestCase):
    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_conductivity(self, mock_opens, mock_json):
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

    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_seebeck(self, mock_opens, mock_json):
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

    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_electronic_thermal_conductivity(self, mock_opens, mock_json):
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

    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_fermi_levels(self, mock_opens, mock_json):
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

    @patch.object(json, 'load')
    @patch("builtins.open", new_callable=mock_open)
    def test_mobility(self, mock_opens, mock_json):
        q = 'mobility'
        data = {q:              {'test': np.zeros((2, 2, 3, 3))},
                'temperatures': [0, 1],
                'doping':       [1, 2]}
        data[q]['test'][0,1,0,0] = 1
        mock_json.return_value = data

        data2 = load.amset('mock', q)
        mock_opens.assert_called_once()
        mock_json.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)
        self.assertEqual(np.array(data2[q])[0,1,0,0,0], 1)
        self.assertEqual(list(data2['scattering_labels']), ['test'])

class AmsetMeshTest(unittest.TestCase):
    @patch.object(h5py, 'File')
    def test_scattering(self, mock_h5py):
        data = {'scattering_rates_up':   np.zeros((1, 2, 2, 3, 1)),
                'scattering_rates_down': np.zeros((1, 2, 2, 3, 1)),
                'scattering_labels':     ['test'],
                'temperatures':        [0, 1],
                'doping':              [1, 2]}
        data['scattering_rates_up'][0,0,1,0,0] = 1
        data['scattering_rates_down'][0,0,1,0,0] = 3
        mock_h5py.return_value = data

        data2 = load.amset_mesh('mock', 'scattering_rates', spin='up')
        mock_h5py.assert_called_once()
        for q2 in ['scattering_rates', 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)
        self.assertEqual(np.array(data2['scattering_rates'])[0,1,0,0,0], 1)

    @patch.object(h5py, 'File')
    def test_spin(self, mock_h5py):
        data = {'scattering_rates_up':   np.zeros((1, 2, 2, 3, 1)),
                'scattering_rates_down': np.zeros((1, 2, 2, 3, 1)),
                'scattering_labels':     ['test'],
                'temperatures':        [0, 1],
                'doping':              [1, 2]}
        data['scattering_rates_up'][0,0,1,0,0] = 1
        data['scattering_rates_down'][0,0,1,0,0] = 3
        mock_h5py.return_value = data

        data2 = load.amset_mesh('mock', 'scattering_rates', spin='avg')
        mock_h5py.assert_called_once()
        for q2 in ['scattering_rates', 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)
        self.assertEqual(np.array(data2['scattering_rates'])[0,1,0,0,0], 2)

class BoltzTraPTest(unittest.TestCase):
    @patch.object(h5py, 'File')
    def test_conductivity(self, mock_h5py):
        q = 'conductivity'
        data = {q:             {'n': np.zeros((2, 2, 3, 3))},
                'temperature': [0, 1],
                'doping':      [1, 2]}
        mock_h5py.return_value = data

        data2 = load.boltztrap('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)

    @patch.object(h5py, 'File')
    def test_effective_mass(self, mock_h5py):
        q = 'average_eff_mass'
        data = {q:             {'n': np.zeros((2, 2, 3, 3))},
                'temperature': [0, 1],
                'doping':      [1, 2]}
        mock_h5py.return_value = data

        data2 = load.boltztrap('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)

    @patch.object(h5py, 'File')
    def test_fermi_level(self, mock_h5py):
        q = 'fermi_level'
        data = {q:             {'n': np.zeros((2, 2, 3, 3))},
                'temperature': [0, 1],
                'doping':      [1, 2]}
        mock_h5py.return_value = data

        data2 = load.boltztrap('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)

    @patch.object(h5py, 'File')
    def test_power_factor(self, mock_h5py):
        q = 'power_factor'
        data = {q:             {'n': np.zeros((2, 2, 3, 3))},
                'temperature': [0, 1],
                'doping':      [1, 2]}
        mock_h5py.return_value = data

        data2 = load.boltztrap('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)

    @patch.object(h5py, 'File')
    def test_seebeck(self, mock_h5py):
        q = 'seebeck'
        data = {q:             {'n': np.zeros((2, 2, 3, 3))},
                'temperature': [0, 1],
                'doping':      [1, 2]}
        mock_h5py.return_value = data

        data2 = load.boltztrap('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)

    @patch.object(h5py, 'File')
    def test_thermal_conductivity(self, mock_h5py):
        q = 'electronic_thermal_conductivity'
        data = {'electronic_thermal_conductivity': {'n': np.zeros((2, 2, 3, 3))},
                'temperature':                     [0, 1],
                'doping':                          [1, 2]}
        mock_h5py.return_value = data

        data2 = load.boltztrap('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'doping', 'meta']:
            self.assertIn(q2, data2)

class Phono3pyTest(unittest.TestCase):
    @patch.object(h5py, 'File')
    def test_frequency(self, mock_h5py):
        q = 'frequency'
        data = {q:             np.zeros((2, 3)),
                'temperature': np.array([0, 1]),
                'qpoint':      np.array([0, 1])}
        mock_data = MagicMock()
        mock_data.__getitem__.side_effect = data.__getitem__
        mock_data.__iter__.side_effect = data.__iter__
        mock_data.__contains__.side_effect = data.__contains__
        mock_h5py.return_value = mock_data

        data2 = load.phono3py('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'meta']:
            self.assertIn(q2, data2)
        mock_data.close.assert_called_once()

    @patch.object(h5py, 'File')
    def test_gamma(self, mock_h5py):
        q = 'gamma'
        data = {q:             np.zeros((2, 2, 3)),
                'temperature': np.array([0, 1]),
                'qpoint':      np.array([0, 1])}
        mock_data = MagicMock()
        mock_data.__getitem__.side_effect = data.__getitem__
        mock_data.__iter__.side_effect = data.__iter__
        mock_data.__contains__.side_effect = data.__contains__
        mock_h5py.return_value = mock_data

        data2 = load.phono3py('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'meta']:
            self.assertIn(q2, data2)
        mock_data.close.assert_called_once()

    @patch.object(h5py, 'File')
    def test_group_velocity(self, mock_h5py):
        q = 'group_velocity'
        data = {q:             np.zeros((2, 3, 3)),
                'temperature': np.array([0, 1]),
                'qpoint':      np.array([0, 1])}
        mock_data = MagicMock()
        mock_data.__getitem__.side_effect = data.__getitem__
        mock_data.__iter__.side_effect = data.__iter__
        mock_data.__contains__.side_effect = data.__contains__
        mock_h5py.return_value = mock_data

        data2 = load.phono3py('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'meta']:
            self.assertIn(q2, data2)
        mock_data.close.assert_called_once()

    @patch.object(h5py, 'File')
    def test_gv_by_gv(self, mock_h5py):
        q = 'gv_by_gv'
        data = {q:             np.zeros((2, 3, 6)),
                'temperature': np.array([0, 1]),
                'qpoint':      np.array([0, 1])}
        mock_data = MagicMock()
        mock_data.__getitem__.side_effect = data.__getitem__
        mock_data.__iter__.side_effect = data.__iter__
        mock_data.__contains__.side_effect = data.__contains__
        mock_h5py.return_value = mock_data

        data2 = load.phono3py('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'meta']:
            self.assertIn(q2, data2)
        mock_data.close.assert_called_once()

    @patch.object(h5py, 'File')
    def test_heat_capacity(self, mock_h5py):
        q = 'heat_capacity'
        data = {q:             np.zeros((2, 2, 3)),
                'temperature': np.array([0, 1]),
                'qpoint':      np.array([0, 1])}
        mock_data = MagicMock()
        mock_data.__getitem__.side_effect = data.__getitem__
        mock_data.__iter__.side_effect = data.__iter__
        mock_data.__contains__.side_effect = data.__contains__
        mock_h5py.return_value = mock_data

        data2 = load.phono3py('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'meta']:
            self.assertIn(q2, data2)
        mock_data.close.assert_called_once()

    @patch.object(h5py, 'File')
    def test_kappa(self, mock_h5py):
        q = 'kappa'
        data = {q:             np.zeros((2, 3)),
                'temperature': np.array([0, 1]),
                'qpoint':      np.array([0, 1])}
        mock_data = MagicMock()
        mock_data.__getitem__.side_effect = data.__getitem__
        mock_data.__iter__.side_effect = data.__iter__
        mock_data.__contains__.side_effect = data.__contains__
        mock_h5py.return_value = mock_data

        data2 = load.phono3py('mock', q)
        mock_h5py.assert_called_once()
        for q2 in ['lattice_thermal_conductivity', 'temperature', 'meta']:
            self.assertIn(q2, data2)
        mock_data.close.assert_called_once()

    @patch.object(h5py, 'File')
    def test_mode_kappa(self, mock_h5py):
        q = 'mode_kappa'
        data = {q:             np.ones((2, 2, 3, 3)),
                'kappa':       np.full((2, 3), 6),
                'weight':      np.ones(2),
                'mesh':        np.full(3, 2),
                'temperature': np.array([0, 1]),
                'qpoint':      np.array([0, 1])}
        mock_data = MagicMock()
        mock_data.__getitem__.side_effect = data.__getitem__
        mock_data.__iter__.side_effect = data.__iter__
        mock_data.__contains__.side_effect = data.__contains__
        mock_h5py.return_value = mock_data

        data2 = load.phono3py('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'meta']:
            self.assertIn(q2, data2)
        mock_data.close.assert_called_once()

    @patch.object(h5py, 'File')
    def test_old_mode_kappa(self, mock_h5py):
        # in old phono3py versions, mode_kappa is multiplied by the mesh
        q = 'mode_kappa'
        data = {q:             np.ones((2, 2, 3, 3)),
                'kappa':       np.full((2, 3), 3),
                'weight':      np.ones(2),
                'mesh':        np.array([2, 1, 1]),
                'temperature': np.array([0, 1]),
                'qpoint':      np.array([0, 1])}
        mock_data = MagicMock()
        mock_data.__getitem__.side_effect = data.__getitem__
        mock_data.__iter__.side_effect = data.__iter__
        mock_data.__contains__.side_effect = data.__contains__
        mock_h5py.return_value = mock_data

        data2 = load.phono3py('mock', q)
        mock_h5py.assert_called_once()
        for q2 in [q, 'temperature', 'meta']:
            self.assertIn(q2, data2)
        mock_data.close.assert_called_once()

    @patch.object(h5py, 'File')
    def test_wideband(self, mock_h5py):
        q = 'wideband'
        data = {'frequency':   np.zeros((2, 3)),
                'gamma':       np.zeros((2, 2, 3)),
                'temperature': np.array([0, 1]),
                'qpoint':      np.array([0, 1])}
        mock_data = MagicMock()
        mock_data.__getitem__.side_effect = data.__getitem__
        mock_data.__iter__.side_effect = data.__iter__
        mock_data.__contains__.side_effect = data.__contains__
        mock_h5py.return_value = mock_data

        data2 = load.phono3py('mock', q)
        mock_h5py.assert_called_once()
        for q2 in ['frequency', 'gamma', 'temperature', 'qpoint', 'meta']:
            self.assertIn(q2, data2)
        mock_data.close.assert_called_once()

class PhonopyDispersionTest(unittest.TestCase):
    @patch.object(yaml, 'safe_load')
    @patch("builtins.open", new_callable=mock_open)
    def test_default(self, mock_opens, mock_yaml):
        mock_data = MagicMock()
        mock_yaml.return_value = mock_data

        data2 = load.phonopy_dispersion('mock')
        mock_opens.assert_called_once()
        mock_yaml.assert_called_once()
        self.assertIn(call('distance'),
                      mock_data['phonon'][0].__getitem__.call_args_list)

    @patch.object(yaml, 'safe_load')
    @patch("builtins.open", new_callable=mock_open)
    def test_scale(self, mock_opens, mock_yaml):
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
    @patch.object(Poscar, 'from_file')
    @patch.object(np, 'loadtxt')
    def test_default(self, mock_load, mock_poscar):
        mock_data = MagicMock()
        mock_load.return_value = mock_data
        mock_POSCAR = MagicMock()
        mock_poscar.return_value = mock_POSCAR

        data2 = load.phonopy_dos('mock')
        mock_load.assert_called_once()
        mock_poscar.assert_called_once()
