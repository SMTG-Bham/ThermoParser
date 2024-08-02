"""Utilities to save data."""

#Functions
#---------
#
#    phono3py:
#        save calculated properties to hdf5.
#    zt:
#        save zt to hdf5 and highlights to yaml.
#    kappa_target:
#        save kappa to hdf5.
#    cumkappa:
#        save cumkappa to text.
#
#    hdf5:
#        save nested dictionaries to hdf5 (up to depth 3).
#    prompt:
#        prompt before overwriting input
#"""

import h5py
import numpy as np
from scipy.interpolate import interp2d
from sys import exit
import tp
import yaml

def shengbte(prefix, output='tp-shengbte', structure='CONTROL',
             atoms=None, kappa_mode='converged'):
    """Saves ShengBTE outputs to a single hdf5 file.

    Also saves dependant properties and metadata.

    Arguments
    ---------

        prefix : str
            ShengBTE data prefix. May include filepath.
        output : str, optional
            output file name (no extension). Default: tp-shengbte.
        structure : str, optional
            structure source. Accepts ShengBTE CONTROL file or VASP POSCAR.
            Overwritten by atoms. Default: CONTROL.
        atoms : str or array-like, optional
            atoms in CONTROL order. Must match number in CONTROL. Overwrites
            structure.
        kappa_mode : str, optional
            method of kappa calculation:

                converged (default):
                    full BTE result.
                rta:
                    relaxation time approximation.
                sg:
                    small grain limit.

    Returns
    -------

        none
            instead writes to hdf5.
    """

    def get_meta(data, q):
        data['meta']['units'][q] = tp.settings.units()[q]
        data['meta']['dimensions'][q] = tp.settings.dimensions()[q]
        return

    ltc = 'lattice_thermal_conductivity'
    kappa_name = {'converged': 'CONV',
                  'rta':       'RTA',
                  'sg':        'sg'}
    data = {'meta': {'kappa_source': 'shengbte_{}'.format(kappa_mode),
                     'units':        {},
                     'dimensions':   {}}}
    
    kdata = np.loadtxt('{}.KappaTensorVsT_{}'.format(prefix, kappa_name[kappa_mode]))
    try:
        data['temperature'] = kdata[:,0].tolist()
        # should probably change kappa to ltc in a future release
        data['kappa'] = np.reshape(kdata[:,1:10], (-1,3,3)).tolist()
    except IndexError:
        data['temperature'] = [kdata[0].tolist()]
        # see above
        data['kappa'] = [np.reshape(kdata[1:10], (3,3)).tolist()]
    get_meta(data, 'temperature')
    get_meta(data, ltc)

    qdata = np.loadtxt('{}.qpoints'.format(prefix))
    data['weight'] = qdata[:,2].tolist()
    get_meta(data, 'weight')
    data['qpoint'] = qdata[:,3:].tolist()
    get_meta(data, 'qpoint')

    data['frequency'] = (np.loadtxt('{}.omega'.format(prefix)) / (2*np.pi)).tolist()
    get_meta(data, 'frequency')
    nq, nbands = np.shape(data['frequency'])

    if atoms is None:
        from pymatgen.io.vasp.inputs import Poscar
        try:
            poscar = Poscar.from_file(structure, check_for_POTCAR=False,
                                      read_velocities=False).as_dict()
            atoms = [p['label'] for p in poscar['structure']['sites']]
        except (FileNotFoundError, ValueError): # CONTROL file
            import re
            with open(structure, 'r') as f:
                atoms = re.search(' elements.*\n', f.read()).group()[10:]
                atoms = atoms.split()
                atoms = [a.strip('"') for a in atoms]
    elif isinstance(atoms, str):
        atoms = atoms.split()
    
    data['dos'] = {}
    data['dos']['total'] = np.loadtxt('{}.dos'.format(prefix))[:,1].tolist()
    pdos = np.loadtxt('{}.pdos'.format(prefix))
    for i, a in enumerate(atoms):
        data['dos'][a] = pdos[:,i+1].tolist()
    get_meta(data, 'dos')

    data['gruneisen'] = np.loadtxt('{}.gruneisen'.format(prefix)).tolist()
    get_meta(data, 'gruneisen')

    v = np.loadtxt('{}.v'.format(prefix))
    v = np.swapaxes(np.reshape(v, (nbands, nq, 3)), 0, 1)
    data['group_velocity'] = v.tolist()
    get_meta(data, 'group_velocity')

    tau = np.loadtxt('{}.w_isotopic'.format(prefix))[:,1]
    tau = np.swapaxes(np.reshape(tau, (nbands, nq)), 0, 1)
    data['lifetime'] = tau.tolist()
    data['gamma'] = (0.5 * np.reciprocal(tau)).tolist()
    get_meta(data, 'lifetime')
    get_meta(data, 'gamma')

    try:
        data['heat_capacity'] = np.loadtxt('{}.cvVsT'.format(prefix))[:,1].tolist()
    except IndexError:
        data['heat_capacity'] = [np.loadtxt('{}.cvVsT'.format(prefix))[1].tolist()]
    get_meta(data, 'heat_capacity')

    hdf5(data, '{}.hdf5'.format(output))

    return

def phono3py(filename, quantities, output='tp-phono3py', force=False):
    """Save calculated properties to hdf5.

    Also saves dependent properties (temperature etc.) and metadata.

    Arguments
    ---------

        filename : str
            filepath.
        quantities : str or list
            values to save. Accepts any phonop3y properties, but only
            lifetime, mean_free_path and/ or occupation are recommended.

        output : str, optional
            output filename (no extension). Default: tp-phono3py.
        force : bool, optional
            force overwrite input file. Default: False.

    Returns
    -------

        none
            instead writes to hdf5.
    """

    if not force:
        prompt(filename, '{}.hdf5'.format(output))
    hdf5(tp.data.load.phono3py(filename, quantities), '{}.hdf5'.format(output))

    return

def zt(efile, kfile=None, direction='avg', doping='n', tinterp=None,
       dinterp=None, kind='linear', output='tp-zt', force=False):
    """Save ZT to hdf5 and highlights to yaml.

    Also saves temperature and doping and metadata.

    Arguments
    ---------

        efile : str
            electronic data filepath.

        kfile : str or float or int
            phononic data filepath or lattice thermal conductivity value.
        direction : str, optional
            crystal direction, accepts x-z/ a-c or average/ avg.
            Default: average.
        doping : str, optional
            doping type for BoltzTraP. Must be n or p. Default: n.

        tinterp : int, optional
            density of interpolation for temperature. None turns it off.
            Default: 200.
        dinterp : int, optional
            density of interpolation for doping. None turns it off.
            Default: 200.
        kind : str, optional
            interpolation kind. Default: linear.

        output : str, optional
            output filename (no extension). Default: tp-zt.
        force : bool, optional
            force overwrite input file. Default: False.

    Returns
    -------

        none
            instead writes to hdf5, yaml and stdout.
    """

    if not force:
        for f in [efile, kfile]:
            prompt(f, ['{}.hdf5'.format(output), '{}.yaml'.format(output)])

    try:
        edata = tp.data.load.amset(efile)
    except UnicodeDecodeError:
        edata = tp.data.load.boltztrap(efile, doping=doping)

    equants = ['conductivity', 'seebeck', 'electronic_thermal_conductivity']
    ltc = 'lattice_thermal_conductivity'
    edata = tp.data.utilities.resolve(edata, equants, direction=direction)

    if kfile is None:
        kfile = 1.
    if not isinstance(kfile, (int, float)):
        kdata = tp.data.load.phono3py(kfile)
        kdata = tp.data.utilities.resolve(kdata, ltc, direction=direction)
        edata, kdata = tp.calculate.interpolate(edata, kdata, 'temperature',
                                                equants, ltc, kind='cubic')
        edata[ltc] = kdata[ltc]
        edata['meta']['dimensions'][ltc] = kdata['meta']['dimensions'][ltc]
        edata['meta']['units'][ltc] = kdata['meta']['units'][ltc]
        edata['meta']['kappa_source'] = kdata['meta']['kappa_source']
    else:
        edata[ltc] = kfile * np.ones(len(edata['temperature']))
        edata['meta']['dimensions'][ltc] = ['temperature']
        edata['meta']['units'][ltc] = tp.settings.units()[ltc]
        edata['meta']['kappa_source'] = 'Set to {} {}'.format(kfile,
                                                  edata['meta']['units'][ltc])
    edata = tp.calculate.zt_fromdict(edata)

    ztdata = {'meta': {**edata['meta'],
                       'original_temperature': np.array(
                                                edata['temperature']).tolist(),
                       'original_doping':      np.array(
                                                edata['doping']).tolist()}}
    # interpolation of zt (if applicable)
    if tinterp is not None or dinterp is not None:
        if tinterp is None:
            ztdata['temperature'] = np.array(edata['temperature']).tolist()
        else:
            ztdata['temperature'] = np.linspace(edata['temperature'][0],
                                                edata['temperature'][-1],
                                                tinterp).tolist()
        if dinterp is None:
            ztdata['doping'] = np.array(edata['doping']).tolist()
        else:
            ztdata['doping'] = np.geomspace(edata['doping'][0],
                                            edata['doping'][-1],
                                            dinterp).tolist()

        ztinterp = interp2d(edata['doping'], edata['temperature'], edata['zt'],
                            kind=kind)
        ztdata['zt'] = ztinterp(ztdata['doping'],
                                ztdata['temperature']).tolist()
    else:
        ztdata['zt'] = np.array(edata['zt']).tolist()
        ztdata['temperature'] = np.array(edata['temperature']).tolist()
        ztdata['doping'] = np.array(edata['doping']).tolist()

    hdf5(ztdata, '{}.hdf5'.format(output))

    # highlights
    ydata = {'meta': ztdata['meta'],
             'max':  {}}

    # max zt and corresponding temperature and doping
    maxindex = np.where(np.round(ztdata['zt'], 10) \
                     == np.round(np.amax(ztdata['zt']), 10))
    ydata['max']['zt'] = ztdata['zt'][maxindex[0][0]][maxindex[1][0]]
    ydata['max']['temperature'] = ztdata['temperature'][maxindex[0][0]]
    ydata['max']['doping'] = ztdata['doping'][maxindex[1][0]]

    # max zt per temperature and corresponding doping
    maxindices = [(i, np.where(np.round(zt, 10) \
                            == np.round(np.amax(zt),10))[0][0]) \
                  for i, zt in enumerate(ztdata['zt'])]
    ydata['zt'] = [ztdata['zt'][i][j] for i, j in maxindices]
    ydata['temperature'] = ztdata['temperature']
    ydata['doping'] = [ztdata['doping'][i] for _, i in maxindices]

    with open('{}.yaml'.format(output), 'w') as f:
        yaml.dump(ydata, f, default_flow_style=False)

    print('Max ZT in the {} direction of {:.2f} at {:.0f} K, {:.2e} carriers cm^-3'.format(
                                                   direction,
                                                   ydata['max']['zt'],
                                                   ydata['max']['temperature'],
                                                   ydata['max']['doping']))

    return

def kappa_target(filename, zt=2, direction='avg', doping='n', tinterp=None,
                 dinterp=None, kind='linear', output='tp-kappa-target',
                 force=False):
    """Save target kappa_l to hdf5.

    Also saves temperature and doping and metadata.

    Arguments
    ---------

        filename : str
            data filepath.

        zt : float, optional
            target ZT. Default: 2.
        direction : str, optional
            crystal direction, accepts x-z/ a-c or average/ avg.
            Default: average.
        doping : str, optional
            doping type for BoltzTraP. Must be n or p. Default: n.

        tinterp : int, optional
            density of interpolation for temperature. None turns it off.
            Default: 200.
        dinterp : int, optional
            density of interpolation for doping. None turns it off.
            Default: 200.
        kind : str, optional
            interpolation kind. Default: linear.

        output : str, optional
            output filename (no extension). Default: tp-kappa-target.
        force : bool, optional
            force overwrite input file. Default: False.

    Returns
    -------

        none
            instead writes to hdf5.
    """

    if not force:
        prompt(filename, '{}.hdf5'.format(output))

    try:
        data = tp.data.load.amset(filename)
    except UnicodeDecodeError:
        data = tp.data.load.boltztrap(filename, doping=doping)
    data = tp.data.utilities.resolve(data, ['conductivity', 'seebeck',
                                   'electronic_thermal_conductivity'],
                                   direction=direction)
    data['zt'] = zt

    data = tp.calculate.kl_fromdict(data)

    kdata = {'meta': {**data['meta'],
                     'original_temperature': data['temperature'],
                     'original_doping':      data['doping']}}
    # interpolation of kl (if applicable)
    if tinterp is not None or dinterp is not None:
        if tinterp is None:
            kdata['temperature'] = data['temperature']
        else:
            kdata['temperature'] = np.linspace(data['temperature'][0],
                                               data['temperature'][-1],
                                               tinterp)
        if dinterp is None:
            kdata['doping'] = data['doping']
        else:
            kdata['doping'] = np.geomspace(data['doping'][0],
                                           data['doping'][-1], dinterp)

        kinterp = interp2d(data['doping'], data['temperature'],
                           data['electronic_thermal_conductivity'], kind=kind)
        kdata['lattice_thermal_conductivity'] = kinterp(kdata['doping'],
                                                        kdata['temperature'])
    else:
        kdata['lattice_thermal_conductivity'] = \
                                       data['lattice_thermal_conductivity']
        kdata['temperature'] = data['temperature']
        kdata['doping'] = data['doping']

    hdf5(kdata, '{}.hdf5'.format(output))

    return

def cumkappa(filename, mfp=False, temperature=300, direction='avg',
             output='tp-cumkappa', extension='dat', force=False):
    """Saves cumulated lattice thermal conductivity against frequency or mfp.

    Saves in normal units and percent to a dat file.

    Arguments
    ---------

        filename : str
            data filepath

        mfp : bool, optional
            calculate against mean free path not frequency. Default: False.
        temperature : float or int, optional
            temperature in K. Default: 300.
        direction : str, optional
            crystal direction, accepts x-z/ a-c or average/ avg.
            Default: average.

        output : str, optional
            output filename (no extension). Default: tp-kappa-target.
        extension : str or list, optional
            output filetype. Must be dat and/ or csv. Default: dat.
        force : bool, optional
            force overwrite input file. Default: False.

    Returns
    -------

        none
            instead writes to dat.
    """

    if isinstance(extension, str):
        extension = [extension]
    csv, dat = False, False
    for e in extension:
        if e == 'csv':
            if not csv:
                csv = True
            else:
                print('Ignoring duplicate extension csv.')
        elif e == 'dat':
            if not dat:
                dat = True
            else:
                print('Ignoring duplicate extension dat.')
        else:
            raise Exception('Extension must be dat and/ or csv.')

    quantity = 'mean_free_path' if mfp else 'frequency'
    data = tp.data.load.phono3py(filename, ['mode_kappa', quantity])
    data = tp.data.utilities.resolve(data, ['mode_kappa', quantity],
                                   temperature=temperature, direction=direction)
    k = np.ravel(data['mode_kappa'])
    q = np.ravel(data[quantity])
    q, k = tp.calculate.cumulate(q, k)
    p = k / np.nanmax(k) * 100

    units = tp.settings.units()
    header = '{}({}) cum_kappa_{d}({}) cum_kappa_{d}(%)'.format(quantity,
                             units[quantity], units['mode_kappa'], d=direction)
    if csv:
        np.savetxt('{}.csv'.format(output), np.transpose([q, k, p]),
                   header=header, delimiter=',')
    if dat:
        np.savetxt('{}.dat'.format(output), np.transpose([q, k, p]),
                   header=header, delimiter=' ')

    return

def hdf5(data, output):
    """Saves to hdf5.

    Aims to make saving nested dictionaries easy, works for 3 layers.

    Arguments
    ---------

        data : dict
            data to save.
        output : str
            output filename.

    Returns
    -------

        None
            instead writes to file.
    """

    with h5py.File(output, 'w') as f:
        for key in data.keys():
            if isinstance(data[key], dict):
                group = f.create_group(key)
                for k in data[key].keys():
                    if isinstance(data[key][k], dict):
                        group2 = group.create_group(k)
                        for k2 in data[key][k].keys():
                            if isinstance(data[key][k][k2], list) and \
                               len(data[key][k][k2]) != 0 and \
                               isinstance(data[key][k][k2][0], str):
                                ds = group2.create_dataset(k2, (len(data[key][k][k2]),),
                                                           dtype=h5py.string_dtype())
                                ds = data[key][k][k2]
                            else:
                                group2[k2] = data[key][k][k2]
                    else:
                        group[k] = data[key][k]
            else:
                f.create_dataset(key, np.shape(data[key]), data=data[key])

    return

def prompt(filename, output):
    """Prompts before overwrite.

    Arguments
    ---------

        filename : str
            input filename.
        output : str or list
            output filename(s).

    Returns
    -------

        none
    """

    if isinstance(output, str):
        output = [output]
    if filename in output:
        tries = 3
        while True:
            print('Warning: this will overwrite {}. Continue?'.format(filename))
            cont = input('[y/n]').lower()
            if cont in ['y', 'ye', 'yes']:
                print('Continuing...')
                break
            elif cont in ['n', 'no']:
                print('Aborting!')
                exit()
            else:
                tries -= 1
                if tries == 2:
                    print('Invalid response. 2 tries remain.')
                elif tries == 1:
                    print('Invalid response. 1 try remains.')
                else:
                    print('Invalid response. Aborting!')
                    exit()

    return
