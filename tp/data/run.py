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
        bta = BoltztrapAnalyzer.from_files(btr_dir)

        norta = {#'average_eff_mass':           bta.get_average_eff_mass,
                 'seebeck':                    bta.get_seebeck}
        rta =   {'conductivity':               bta.get_conductivity,}
                 #'power_factor':               bta.get_power_factor,
                 #'thermal_conductivity':       bta.get_thermal_conductivity}
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
                'temperature':                temperature,}
                #'thermal_conductivity':       {}}

        for d in ['n', 'p']:
            for key in norta.keys():
                data[key][d] = [norta[key](output='tensor')[d][temperature[i]] \
                    for i in range(len(temperature))]
            for key in rta.keys():
                data[key][d] = [rta[key](output='tensor',
                    relaxation_time=relaxation_time)[d][temperature[i]] \
                    for i in range(len(temperature))]
            #for key in lmbda.keys():
            #    data[key][d] = [lmbda[key](output='tensor',
            #        temp=temperature[i], Lambda=Lambda, doping_levels=True)[d] \
            #        for i in range(len(temperature))]
        data['mu_bounds'] = [bta.get_mu_bounds(temp=temperature[i]) \
            for i in range(len(temperature))]
        data['hall_carrier_concentration'] = \
            [bta.get_hall_carrier_concentration()[temperature[i]] \
            for i in range(len(temperature))]

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
