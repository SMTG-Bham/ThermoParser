"""External code running tools."""

#Functions
#---------
#
#    boltztrap
#    boltztrap2
#"""

import numpy as np
import shutil
import tp

workers = tp.settings.get_workers()

def boltztrap(tmax=1001, tstep=50, tmin=None, doping=np.logspace(18, 21, 17),
              ke_mode='boltzmann', vasprun='vasprun.xml', kpoints=None,
              relaxation_time=1e-14, lpfac=10, run=True, analyse=True,
              output='boltztrap.hdf5', run_dir='.', clean=False, **kwargs):
    """Runs BoltzTraP from a VASP density of states (DoS).

    Wrapper for pymatgen.electronic_structure.boltztrap but runs faster than 
    using the built in from_files method and outputs an hdf5 file.
    Note: BoltzTraP can be a fickle friend, so if you're getting errors,
    it may be worth reinstalling or trying on a different machine.

    Arguments
    ---------

        tmax : float, optional
            maximum temperature in K. Default: 1000.
        tstep : float, optional
            temperature step in K. Default: 50.
        tmin : float, optional
            minimum temperature in K. This does not reduce how many
            temperatures are run in BoltzTraP, only how many are saved
            to hdf5. Default: tstep.
        doping : array-like, optional
            doping concentrations in cm-1.
            Default: np.logspace(18, 21, 17).
        ke_mode : str, optional
            method for calculating the electronic thermal conductivity.
            Options:

                boltzmann (default):
                    standard boltztrap method. Madsen and Singh,
                    Comput. Phys. Commun. 2006, 175, 67.
                wiedemann:
                    Wiedemann-Franz law with constant L = 2.44E-8.
                    Franz and Wiedemann, Ann. Phys. 1853, 165, 497.
                snyder:
                    Wiedemann-Franz law, with L varying with Seebeck.
                    Kim et al., APL Mat. 2015, 3, 041506.

        vasprun : str, optional
            path to vasprun. Default: vasprun.xml.
        kpoints : str, optional
            path to KPOINTS file if there are zero-weighted k-points.
            Default: KPOINTS.
        relaxation_time : float, optional
            charge carrier relaxation time. Default: 1e-14.
        lpfac : int, optional
            DoS interpolation factor. Default: 10.

        run : bool, optional
            run BoltzTraP. Default: True.
        analyse : bool, optional
            analyse BoltzTraP. Default: True.
        output : str, optional
            output hdf5 filename. Default: boltztrap.hdf5.
        run_dir : str, optional
            path to run boltztrap in. Default: current directory.
        clean : bool, optional
            remove boltztrap directory post-run. Default: False.

        kwargs
            passed to pymatgen.electronic.structure.boltztrap.BoltztrapRunner.

    Returns
    -------

        None
            instead prints to hdf5 (see below).

    hdf5 File Contents
    ------------------

        average_eff_mass : dict
            the charge carrier effective mass in units of m_e, taking
            into account all bands, as opposed to the single parabolic
            band model.
            Data is a dictionary with an array each for n and p doping,
            of shape (temperatures, concentrations, 3, 3).
        conductivity : dict
            electric conductivity in S m-1.
            Data is a dictionary with an array each for n and p doping,
            of shape (temperatures, concentrations, 3, 3).
        doping : array-like
            carrier concentration in cm-1. Identical to input.
        electronic_thermal_conductivity : dict
            electronic thermal conductivity in W m-1 K-1.
            Data is a dictionary with an array each for n and p doping,
            of shape (temperatures, concentrations, 3, 3).
        fermi_level : dict
            fermi level at different temperatures in units of eV.
            Data is a dictionary with an array each for n and p doping,
            of shape (temperatures, concentrations).
        power_factor : dict
            power factor in W m-1 K-2.
            Data is a dictionary with an array each for n and p doping,
            of shape (temperatures, concentrations, 3, 3).
        seebeck : dict
            Seebeck coefficient in muV K-1.
            Data is a dictionary with an array each for n and p doping,
            of shape (temperatures, concentrations, 3, 3).
        temperature : numpy array
            temperatures in K.
        meta : dict
            metadata:

                dimensions : dict
                    dimensions of each property above.
                interpolation_factor : int
                    lpfac.
                ke_mode : str
                    as input.
                relaxation_time : float
                    as input.
                soc : bool
                    spin-orbit coupling calculation.
                units : dict
                    units of each property above.
    """

    import os
    from pymatgen.electronic_structure.boltztrap \
         import BoltztrapRunner, BoltztrapAnalyzer, BoltztrapError
    from pymatgen.io.vasp.outputs import Vasprun
    from scipy import constants

    # check inputs

    for name, value in zip(['run', 'analyse', 'clean'],
                           [ run,   analyse,   clean]):
        assert isinstance(value, bool), '{} must be True or False'.format(name)

    ke_mode = ke_mode.lower()
    ke_modes =  ['boltzmann', 'wiedemann', 'snyder']
    assert ke_mode in ke_modes, 'ke_mode must be {} or {}.'.format(
                                ', '.join(ke_modes[:-1]), ke_modes[-1])

    run_dir = os.path.abspath(run_dir)

    tmax += tstep
    tmin = tstep if tmin is None else tmin
    temperature = np.arange(tmin, tmax, tstep)

    if run: # run boltztrap from vasprun.xml -> boltztrap directory
        doping = np.array(doping)
        vr = Vasprun(vasprun, parse_potcar_file=False)
        soc = vr.as_dict()['input']['parameters']['LSORBIT']
        try:
            bs = vr.get_band_structure(line_mode=False)
            nelect = vr.parameters['NELECT']
            btr = BoltztrapRunner(bs, nelect, doping=list(doping), tmax=tmax,
                                  tgrid=tstep, soc=soc, lpfac=lpfac, **kwargs)
            print('Running Boltztrap...', end='')
            btr_dir = btr.run(path_dir=run_dir)
            print('Done.')
        except BoltztrapError:
            bs = vr.get_band_structure(line_mode=True,kpoints_filename=kpoints)
            nelect = vr.parameters['NELECT']
            btr = BoltztrapRunner(bs, nelect, doping=list(doping), tmax=tmax,
                                  tgrid=tstep, soc=soc, lpfac=lpfac, **kwargs)
            btr_dir = btr.run(path_dir=run_dir)
            print('Done.')

        """
        Detects whether the BoltzTraP build on this computer writes the
        doping concentrations correctly, and if it doesn't, writes them.
        """
        with open(os.path.join(btr_dir, 'boltztrap.outputtrans'), 'r') as f:
            line = f.readlines()[-2]
        if len(line) >= 23 and line[:23] == ' Calling FermiIntegrals':
            with open(os.path.join(btr_dir, 'boltztrap.outputtrans'),'a') as f:
                for i, x in enumerate(np.concatenate((doping, -doping))):
                    f.write(
                   'Doping level number {} n = {} carriers/cm3\n'.format(i, x))
    else:
        btr_dir = os.path.join(run_dir, 'boltztrap')

    if analyse: # run boltztrap from boltztrap directory -> hdf5 file
        print('Analysing Boltztrap...', end='')
        bta = BoltztrapAnalyzer.from_files(btr_dir)
        print('Done.')

        data = {'average_eff_mass':                {},
                'conductivity':                    {},
                'doping':                          doping,
                'electronic_thermal_conductivity': {},
                'power_factor':                    {},
                'seebeck':                         {},
                'temperature':                     temperature,
                'meta':
                    {'units':      {'average_eff_mass':                'm_e',
                                    'conductivity':                    'S m-1',
                                    'doping':                          'cm-1',
                                    'electronic_thermal_conductivity': 'W m-1 K-1',
                                    'fermi_level':                     'eV',
                                    'power_factor':                    'W m-1 K-2',
                                    'seebeck':                         'muV K-1',
                                    'temperature':                     'K'},
                     'dimensions': {'average_eff_mass':                ['temperature', 'doping', 3, 3],
                                    'conductivity':                    ['temperature', 'doping', 3, 3],
                                    'doping':                          ['doping'],
                                    'electronic_thermal_conductivity': ['temperature', 'doping', 3, 3],
                                    'fermi_level':                     ['temperature', 'doping'],
                                    'power_factor':                    ['temperature', 'doping', 3, 3],
                                    'seebeck':                         ['temperature', 'doping', 3, 3],
                                    'temperature':                     ['temperature']},
                   'interpolation_factor': lpfac,
                   'ke_mode':              ke_mode,
                   'relaxation_time':      relaxation_time,
                   'soc':                  soc}}

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
            data['power_factor'][d] = tp.calculate.power_factor(
                                                   data['conductivity'][d],
                                                   data['seebeck'][d],
                                                   use_tprc=False)
            if ke_mode == 'boltztrap':
                data['electronic_thermal_conductivity'][d] = \
                      k[d] * relaxation_time \
                      - data['power_factor'][d] * temperature[:,None,None,None]
            else:
                if ke_mode == 'wiedemann':
                    L = 2.44E-8
                else:
                    L = (1.5 + np.exp(-np.abs(data['seebeck'][d])/116)) * 1e-8
                data['electronic_thermal_conductivity'][d] = \
                    L * data['conductivity'][d] * temperature[:,None,None,None]

        data['fermi_level'] = f
        print('Done.')

        tp.data.save.hdf5(data, output)

        if clean:
            shutil.rmtree(btr_dir)
            os.remove(os.path.join(run_dir, 'boltztrap.out'))

    return

@tp.docstring_replace(workers=str(workers))
def boltztrap2(tmax=1001, tstep=50, tmin=None, doping=np.logspace(18, 21, 22),
               scissor=None, relaxation_time=1e-14, interpolate=True, lpfac=10,
               emin=-0.2, emax=0.2, dosbins=4000, run=True, prefix='tp',
               dos_dir='.', run_dir='.', workers=workers):
    """Runs BoltzTraP2 from a density of states (DoS).

    Heavily based on the official tutorial and CLI dope function, but
    outputs to a tp hdf5 file.

    Arguments
    ---------

        tmax : float, optional
            maximum temperature in K. Default: 1000.
        tstep : float, optional
            temperature step in K. Default: 50.
        tmin : float, optional
            minimum temperature in K. Default: tstep.
        doping : array-like, optional
            doping concentrations in cm-1.
            Default: np.logspace(18, 21, 22).
        scissor : float, optionao
            target band gap in eV. Default: None.

        relaxation_time : float, optional
            charge carrier relaxation time. Default: 1e-14.
        interpolate : bool, optional
            run DoS interpolation. Default: True.
        lpfac : int, optional
            DoS interpolation factor. Default: 10.
        emin : float, optional
            minimum energy from the Fermi level to include bands during
            interpolation. Default: -0.2.
        emax : float, optional
            maximum energy from the Fermi level to include bands during
            interpolation. Default: +0.2.
        dosbins : int, optional
            number of energy bins to sort the DoS into. Default: 4000.

        run : bool, optional
            run BoltzTraP2. Default: True.
        prefix : str, optional
            output filename prefix. Default: tp.
        dos_dir : str, optional
            path to dft DoS output. Default: current directory.
        run_dir : str, optional
            path to run boltztrap in. Default: current directory.
        workers : int, optional
            number of workers for paralellised section.
            Default: {workers}.

    Returns
    -------

        None
            instead prints to hdf5 (see below).

    hdf5 File Contents
    ------------------

        conductivity : array-like
            electric conductivity in S m-1.
            Shape: (temperatures, concentrations, 3, 3).
        doping : array-like
            carrier concentration in cm-1. Identical to input.
        electronic_heat_capacity : array-like
            electronic contributions to the specific heat capacity in
            J mol-1 K-1. Shape (temperatures, concentrations)
        electronic_thermal_conductivity : array-like
            electronic thermal conductivity in W m-1 K-1.
            Shape (temperatures, concentrations, 3, 3).
        fermi_level : array-like
            fermi level at different temperatures in eV.
            Shape: (temperatures, concentrations).
        hall : array-like
            Hall coefficient in J C-1.
            Shape: (temperatures, concentrations, 3, 3, 3).
        power_factor : array-like
            power factor in W m-1 K-2.
            Shape: (temperatures, concentrations, 3, 3).
        seebeck : dict
            Seebeck coefficient in muV K-1.
            Shape: (temperatures, concentrations, 3, 3).
        temperature : numpy array
            temperatures in K.
        meta : dict
            metadata:

                dimensions : dict
                    dimensions of each property above.
                interpolation_factor : int
                    lpfac.
                relaxation_time : float
                    as input.
                units : dict
                    units of each property above.
    """

    import os
    from BoltzTraP2 import bandlib as bl
    from BoltzTraP2 import dft
    from BoltzTraP2 import fite, serialization, sphere, units

    # check inputs

    for name, value in zip(['run', 'interpolate'],
                           [ run,   interpolate]):
        assert isinstance(value, bool), '{} must be True or False'.format(name)

    run_dir = os.path.abspath(run_dir)
    dos_dir = os.path.abspath(dos_dir)
    ifile = os.path.join(dos_dir, f'{prefix}.btp2')

    tmax += tstep
    tmin = tstep if tmin is None else tmin
    temperature = np.arange(tmin, tmax, tstep)

    if interpolate:
        print('Loading DFT data...')
        data = dft.DFTData(dos_dir)
        print('Interpolating DFT data...')
        data.bandana(emin=data.fermi+emin, emax=data.fermi+emax)
        eq = sphere.get_equivalences(data.atoms, data.magmom,
                                     nkpt=len(data.kpoints) * lpfac)
        coeffs = fite.fitde3D(data, eq)
        meta = serialization.gen_bt2_metadata(data, data.mommat is not None)
        serialization.save_calculation(ifile, data, eq, coeffs, meta)
    else:
        print('Loading interpolated DoS...')
        data, eq, coeffs, meta = serialization.load_calculation(ifile)

    if run:
        print('Running BoltzTraP2...')
        lattvec = data.get_lattvec()
        eband, vvband, cband = fite.getBTPbands(eq, coeffs, lattvec,
                                                curvature=True, nworkers=workers)
        epsilon, dos, vvdos, cdos = bl.BTPDOS(eband, vvband, cband,
                                              npts=dosbins, Tmin=tmin)
        if scissor is not None:
            scissor *= units.eV
            eband = bl.apply_scissor(epsilon, dos, data.nelect, eband, scissor,
                                     data.dosweight)
            epsilon, dos, vvdos, cdos = bl.BTPDOS(eband, vvband, cband,
                                                  npts=dosbins, tmin=tmin)
        margin = 9. * units.BOLTZMANN * temperature.max()
        mumin = epsilon.min() + margin
        mumax = epsilon.max() - margin

        v = data.atoms.get_volume() * units.Angstrom**3
        vcm = data.atoms.get_volume() * 1e-24
        mu = np.empty((len(temperature), len(doping)))
        #sdos = np.empty((len(temperature), len(doping)))
        cv = np.empty((len(temperature), len(doping)))
        n = np.empty((len(temperature), len(doping)))
        l0 = np.empty((len(temperature), len(doping), 3, 3))
        l1 = np.empty((len(temperature), len(doping), 3, 3))
        l2 = np.empty((len(temperature), len(doping), 3, 3))
        lm11 = np.empty((len(temperature), len(doping), 3, 3, 3))
        sigma = np.empty((len(temperature), len(doping), 3, 3))
        seebeck = np.empty((len(temperature), len(doping), 3, 3))
        kappa = np.empty((len(temperature), len(doping), 3, 3))
        hall = np.empty((len(temperature), len(doping), 3, 3, 3))
        for i, t in enumerate(temperature):
            dmin = (bl.calc_N(epsilon, dos, mumax, t, data.dosweight) \
                 +  data.nelect) / vcm
            dmax = (bl.calc_N(epsilon, dos, mumin, t, data.dosweight) \
                 +  data.nelect) / vcm
            if doping.min() <= dmin or doping.max() >= dmax:
                raise ValueError(f'doping must be between {dmin:.2g} and '
                                 f'{dmax:.2g} at {t:.2f}.')

            for j, d in enumerate(doping):
                ne = data.nelect - d * vcm
                mu[i, j] = bl.solve_for_mu(epsilon, dos, ne, t, data.dosweight,
                                           refine=True, try_center=False)
                
            #sdos[i,:] = bl.smoothen_DOS_direct(epsilon, dos, t, mu[i, :])
            cv[i] = bl.calc_cv(epsilon, dos, mu[i,:], np.array([t]),
                               dosweight=data.dosweight)
            n[i], l0[i], l1[i], l2[i], lm11[i] = bl.fermiintegrals(
                epsilon, dos, vvdos, mur=mu[i], Tr=np.array([t]),
                dosweight=data.dosweight, cdos=cdos)
        sigma, seebeck, kappa, hall = bl.calc_Onsager_coefficients(
                                          l0, l1, l2, mu, temperature, v, lm11)
        sigma *= relaxation_time
        kappa *= relaxation_time
        seebeck *= 1e6

        data = {'conductivity':                    sigma,
                'doping':                          doping,
                'electronic_heat_capacity':        cv,
                'electronic_thermal_conductivity': kappa,
                'fermi_level':                     mu,
                'hall':                            hall,
                'power_factor':                    sigma*seebeck**2,
                'seebeck':                         seebeck,
                'temperature':                     temperature,
                'meta':
                    {'units':      {'conductivity':                    'S m-1',
                                    'doping':                          'cm-1',
                                    'electronic_heat_capacity':        'J mol-1 K-1',
                                    'electronic_thermal_conductivity': 'W m-1 K-1',
                                    'fermi_level':                     'eV',
                                    'hall':                            'm3 C-1',
                                    'power_factor':                    'W m-1 K-2',
                                    'seebeck':                         'muV K-1',
                                    'temperature':                     'K'},
                     'dimensions': {'conductivity':                    ['temperature', 'doping', 3, 3],
                                    'doping':                          ['doping'],
                                    'electronic_heat_capacity':        ['temperature', 'doping'],
                                    'electronic_thermal_conductivity': ['temperature', 'doping', 3, 3],
                                    'fermi_level':                     ['temperature', 'doping'],
                                    'hall':                            ['temperature', 'doping', 3, 3, 3],
                                    'power_factor':                    ['temperature', 'doping', 3, 3],
                                    'seebeck':                         ['temperature', 'doping', 3, 3],
                                    'temperature':                     ['temperature']},
                   'interpolation_factor': lpfac,
                   'relaxation_time':      relaxation_time}}

        print('Done.')
        tp.data.save.hdf5(data, f'{prefix}-btp2.hdf5')

    return
