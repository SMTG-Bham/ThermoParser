"""Code running tools.

Functions:
    boltztrap
"""

import numpy as np
import tp

def boltztrap(tmax=1000, tstep=10, doping=np.logspace(18, 21, 100),
              vasprun='vasprun.xml', soc=False, zero_weighted=True, lpfac=10,
              relaxation_time=1e-14, Lambda=0.5, run=True, analyse=True,
              output='boltztrap.hdf5'):
    """Runs BoltzTraP

    Also writes to a file.
    Some quantities require the fdint package, which is
    pip-installable.
    Note: BoltzTraP can be a fickle friend, so if you're getting
    errors, it may be worth reinstalling or trying on a different
    machine.

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
            zero weighted kpoints used. Default: True.
        lpfac : int, optional
            interpolate the DoS k-points by lpfac times. Default: 10.
        relaxation_time : float, optional
            charge carrier relaxation time. Default: 1e-14.
        Lambda : float, optional
            scattering fitting parameter. Default: 0.5.

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

    temperature = np.arange(200, tmax, tstep)

    if run:
        doping = np.array(doping)
        vr = Vasprun(vasprun)
        bs = vr.get_band_structure(line_mode=zero_weighted)
        nelect = vr.parameters['NELECT']

        btr = BoltztrapRunner(bs, nelect, doping=list(doping), tmax=tmax, tgrid=tstep,
                              lpfac=lpfac)
        print('Running Boltztrap...', end='')
        btr_dir = btr.run(path_dir='.')
        print('Done.')

        with open(os.path.join(btr_dir, 'boltztrap.outputtrans'), 'a') as f:
            for i, x in enumerate(np.concatenate((doping, -doping))):
                f.write(
                  'Doping level number {} n = {} carriers/cm3\n'.format(i, x))
    else:
        btr_dir = 'boltztrap'

    if analyse:
        from time import perf_counter
        start=perf_counter()
        bta = BoltztrapAnalyzer.from_files(btr_dir)

        data = {#'average_eff_mass':           {},
                #'complexity_factor':          {},
                'conductivity':               {},
                'doping':                     doping,
                'power_factor':               {},
                'seebeck':                    {},
                #'seebeck_eff_mass':           {},
                'temperature':                temperature,
                'thermal_conductivity':       {}}

        c = bta._cond_doping
        s = bta._seebeck_doping
        k = bta._kappa_doping

        for d in ['n', 'p']:
            c[d] = np.array([c[d][t] for t in temperature])
            s[d] = np.array([s[d][t] for t in temperature])
            k[d] = np.array([k[d][t] for t in temperature])
            data['conductivity'][d] = np.multiply(c[d], relaxation_time)
            data['seebeck'][d] = np.multiply(s[d], 1e6)
            data['power_factor'][d] = tp.calculate.power_factor(c[d], s[d])
            data['thermal_conductivity'][d] = (k[d] - s[d] ** 2 * c[d] \
                                            * temperature[:,None,None,None]) \
                                            * relaxation_time

        #sem = [bta.get_seebeck_eff_mass(output='tensor', temp=t, \
        #       Lmbda=Lambda, doping_levels=True) for t in temperature]
        #print('sem', perf_counter() - start)
        #cf = [bta.get_complexity_factor(output='tensor', temp=t, \
        #      Lambda=Lambda, doping_levels=True) for t in temperature]
        #print('cf', perf_counter() - start)
        #data['mu_bounds'] = [bta.get_mu_bounds(temp=temperature[i]) \
        #                     for i in range(len(temperature))]
        #print('mu', perf_counter() - start)
        #data['hall_carrier_concentration'] = \
        #      [bta.get_hall_carrier_concentration()[temperature[i]] \
        #       for i in range(len(temperature))]
        #print('hall', perf_counter() - start)

        #for d in ['n', 'p']:
        #    data['seebeck_eff_mass'][d] = [s[d] for s in sem]
        #    data['complexity_factor'][d] = [c[d] for c in cf]

        datafile = h5py.File(output, 'w')
        for key in data.keys():
            if isinstance(data[key], dict):
                group = datafile.create_group(key)
                for k in data[key].keys():
                    group[k] = data[key][k]
            else:
                datafile.create_dataset(key, np.shape(data[key]), data=data[key])
        datafile.close()

    return

def boltztrap_serial(tmax=1000, tstep=10, doping=np.logspace(18, 21, 100),
              vasprun='vasprun.xml', soc=False, zero_weighted=True, lpfac=10,
              relaxation_time=1e-14, Lambda=0.5, run=True, analyse=True,
              output='boltztrap.hdf5'):
    """Runs BoltzTraP

    Also writes to a file.

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
            zero weighted kpoints used. Default: True.
        lpfac : int, optional
            interpolate the DoS k-points by lpfac times. Default: 10.
        relaxation_time : float, optional
            charge carrier relaxation time. Default: 1e-14.
        Lambda : float, optional
            scattering fitting parameter. Default: 0.5.

        run : bool, optional
            run BoltzTraP. Default: True.
        analyse : bool, optional
            analyse BoltzTraP. Default: True.
        output : str, optional
            output hdf5 filename. Default: boltztrap.hdf5.
    """

    import h5py
    import os
    from pymatgen.electronic_structure.boltztrap import BoltztrapRunner, BoltztrapAnalyzer
    from pymatgen.io.vasp.outputs import Vasprun

    temperature = np.arange(200, tmax, tstep)

    if run:
        doping = np.array(doping)
        vr = Vasprun(vasprun)
        bs = vr.get_band_structure(line_mode=zero_weighted)
        nelect = vr.parameters['NELECT']

        btr = BoltztrapRunner(bs, nelect, doping=list(doping), tmax=tmax, tgrid=tstep,
                              lpfac=lpfac)
        btr_dir = btr.run(path_dir='.')

        with open(os.path.join(btr_dir, 'boltztrap.outputtrans'), 'a') as f:
            for i, x in enumerate(np.concatenate((doping, -doping))):
                f.write(
                  'Doping level number {} n = {} carriers/cm3\n'.format(i, x))
    else:
        btr_dir = 'boltztrap'

    if analyse:
        from time import perf_counter
        start = perf_counter()
        bta = BoltztrapAnalyzer.from_files(btr_dir)
        print('bta', perf_counter() - start)

        norta = {#'average_eff_mass':           bta.get_average_eff_mass,
                 'seebeck':                    bta.get_seebeck}
        rta =   {'conductivity':               bta.get_conductivity,
                 #'power_factor':               bta.get_power_factor,
                 'thermal_conductivity':       bta.get_thermal_conductivity}
        lmbda = {'seebeck_eff_mass':           bta.get_seebeck_eff_mass,
                 'complexity_factor':          bta.get_complexity_factor}
        other = {'mu_bounds':                  bta.get_mu_bounds,
                 'hall_carrier_concentration': bta.get_hall_carrier_concentration}

        data = {#'average_eff_mass':           {},
                #'complexity_factor':          {},
                'conductivity':               {},
                'doping':                     doping,
                #'power_factor':               {},
                'seebeck':                    {},
                #'seebeck_eff_mass':           {},
                'temperature':                temperature,
                'thermal_conductivity':       {}}

        for d in ['n', 'p']:
            for key in norta.keys():
                data[key][d] = [norta[key](output='tensor')[d][temperature[i]] \
                    for i in range(len(temperature))]
                print(key, perf_counter() - start)
            for key in rta.keys():
                data[key][d] = [rta[key](output='tensor',
                    relaxation_time=relaxation_time)[d][temperature[i]] \
                    for i in range(len(temperature))]
                print(key, perf_counter() - start)
            #for key in lmbda.keys():
            #    data[key][d] = [lmbda[key](output='tensor',
            #        temp=temperature[i], Lambda=Lambda, doping_levels=True)[d] \
            #        for i in range(len(temperature))]
            #    print(key, perf_counter() - start)
        #data['mu_bounds'] = [bta.get_mu_bounds(temp=temperature[i]) \
        #    for i in range(len(temperature))]
        #print('mu', perf_counter() - start)
        #data['hall_carrier_concentration'] = \
        #    [bta.get_hall_carrier_concentration()[temperature[i]] \
        #    for i in range(len(temperature))]
        #print('hall', perf_counter() - start)

        datafile = h5py.File(output, 'w')
        for key in data.keys():
            if isinstance(data[key], dict):
                group = datafile.create_group(key)
                for k in data[key].keys():
                    group[k] = data[key][k]
            else:
                datafile.create_dataset(key, np.shape(data[key]), data=data[key])
        datafile.close()
        print('write', perf_counter() - start)

    return
