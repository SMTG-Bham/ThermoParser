"""Code running tools.

Functions:
    boltztrap
"""

import numpy as np
import tp

def boltztrap(tmax=1000, tstep=10, doping=np.logspace(18, 21, 100),
              vasprun='vasprun.xml', soc=False, zero_weighted=False, lpfac=10,
              relaxation_time=1e-14, run=True, analyse=True,
              output='boltztrap.hdf5'):
    """Runs BoltzTraP

    Runs quicker than the pymatgen from_files version.
    Also writes to a hdf5 file.
    Note: BoltzTraP can be a fickle friend, so if you're getting errors,
    it may be worth reinstalling or trying on a different machine.

    Arguments:
        tmax : float, optional
            maximum temperature in K. Default: 1000.
        tstep : float, optional
            temperature step in K. Default: 10.
        doping : array-like, optional
            doping concentrations in cm-1. Default: np.logspace(18, 21, 101).

        vasprun : str, optional
            path to vasprun. Default: vasprun.xml.
        soc : bool, optional
            spin orbit coupling used. Default: False.
        zero_weighted : bool, optional
            zero weighted kpoints used. Default: False.
        lpfac : int, optional
            interpolate the DoS k-points by lpfac times. Default: 10.
        relaxation_time : float, optional
            charge carrier relaxation time. Default: 1e-14.

        run : bool, optional
            run BoltzTraP. Default: True.
        analyse : bool, optional
            analyse BoltzTraP. Default: True.
        output : str, optional
            output hdf5 filename. Default: boltztrap.hdf5.
    """

    import h5py
    import os
    from multiprocessing import Pool
    from pymatgen.electronic_structure.boltztrap import BoltztrapRunner, BoltztrapAnalyzer
    from pymatgen.io.vasp.outputs import Vasprun
    from scipy import constants

    temperature = np.arange(200, tmax, tstep)

    if run: # run boltztrap from vasprun.xml -> boltztrap directory
        doping = np.array(doping)
        vr = Vasprun(vasprun)
        bs = vr.get_band_structure(line_mode=zero_weighted)
        nelect = vr.parameters['NELECT']

        btr = BoltztrapRunner(bs, nelect, doping=list(doping), tmax=tmax,
                                  tgrid=tstep, lpfac=lpfac)
        print('Running Boltztrap...', end='')
        btr_dir = btr.run(path_dir='.')
        print('Done.')

        with open(os.path.join(btr_dir, 'boltztrap.outputtrans'), 'a') as f:
            for i, x in enumerate(np.concatenate((doping, -doping))):
                f.write(
                  'Doping level number {} n = {} carriers/cm3\n'.format(i, x))
    else:
        btr_dir = 'boltztrap'

    if analyse: # run boltztrap from boltztrap directory -> hdf5 file
        print('Analysing Boltztrap...', end='')
        bta = BoltztrapAnalyzer.from_files(btr_dir)
        print('Done.')

        data = {'average_eff_mass':     {},
                'conductivity':         {},
                'doping':               doping,
                'power_factor':         {},
                'seebeck':              {},
                'temperature':          temperature,
                'thermal_conductivity': {}}

        # load data

        print('Calculating transport...', end='')
        c = bta._cond_doping
        dp = bta.doping
        e = constants.e
        f = bta.mu_doping
        k = bta._kappa_doping
        me = constants.m_e
        s = bta._seebeck_doping

        # calculate transport properties

        for d in ['n', 'p']:
            c[d] = np.array([c[d][t] for t in temperature])
            dp[d] = np.array(dp[d])
            k[d] = np.array([k[d][t] for t in temperature])
            f[d] = np.array([f[d][t] for t in temperature])
            s[d] = np.array([s[d][t] for t in temperature])
            data['average_eff_mass'][d] = np.linalg.inv(c[d]) \
                                        * dp[d][None, :, None, None] \
                                        * 1e6 * e ** 2 / me
            data['conductivity'][d] = np.multiply(c[d], relaxation_time)
            data['seebeck'][d] = np.multiply(s[d], 1e6)
            data['power_factor'][d] = tp.calculate.power_factor(c[d], s[d])
            data['thermal_conductivity'][d] = (k[d] - s[d] ** 2 * c[d] \
                                            * temperature[:,None,None,None]) \
                                            * relaxation_time
        data['fermi_level'] = f
        print('Done.')

        # write data

        datafile = h5py.File(output, 'w')
        for key in data.keys():
            if isinstance(data[key], dict):
                group = datafile.create_group(key)
                for k in data[key].keys():
                    group[k] = data[key][k]
            else:
                datafile.create_dataset(key, np.shape(data[key]),
                                        data=data[key])
        datafile.close()

    return
