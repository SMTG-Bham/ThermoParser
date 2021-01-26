"""Data loading tools.

Loads data from codes into a dictionary, with units and array structures
standardised. Also adds a ``meta`` subdictionary, which contains units
and the data source.

Functions
---------

    amset
    amset_mesh
    boltztrap
    phono3py
    phonopy_dispersion
    phonopy_dos

    get_path:
        gets high path from phonopy dispersion data.
"""

import numpy as np
import tp
from tp import settings

def amset(filename, quantities=['temperature', 'doping', 'seebeck',
          'conductivity', 'electronic_thermal_conductivity']):
    """Loads AMSET transport data.

    Includes unit conversion and outputs units (see tp.settings).
    Swaps temperature and doping indices so temperature is first, for
    consistency with other codes, and mobility is formatted as an array,
    like scattering_rates is in the mesh.h5 file. Maintains basic
    compatibility with amset 0.1.

    Arguments
    ---------

        filename : str
            filepath.

        quantites : dict, optional
            values to extract. Default: temperature, doping, seebeck,
            conductivity, electronic_thermal_conductivity.

    Returns
    -------

        dict
            extracted values.
    """

    import json

    # name conversions and abbreviations

    conversions = settings.amset_conversions()
    anames = settings.to_amset()
    tnames = settings.to_tp()
    units = settings.units()
    if isinstance(quantities, str): quantities = quantities.split()
    quantities = [anames[q] if q in anames else q for q in quantities]

    # list of quantities dependant on doping, temperature and scattering
    hasdope = ['fermi_levels', 'conductivity', 'seebeck',
               'electronic_thermal_conductivity', 'mobility']
    hastemp = ['fermi_levels', 'conductivity', 'seebeck',
               'electronic_thermal_conductivity', 'mobility']
    hastype = ['mobility']

    # add dependant variables
    if 'doping' not in quantities:
        for q in quantities:
            if q in hasdope:
                quantities.append('doping')
                break
    if 'temperatures' not in quantities:
        for q in quantities:
            if q in hastemp:
                quantities.append('temperatures')
                break

    # load data

    with open(filename) as f:
        data = json.load(f)

    data2 = {'meta': {'electronic_source': 'amset',
                      'units':             {}}}
    for q in quantities:
        assert q in data, \
               '{} unrecognised. Quantity must be in {} or {}.'.format(q,
               ', '.join(list(data)[:-1]), list(data)[-1])
        q2 = tnames[q] if q in tnames else q
        # compatibility with previous version
        if isinstance(data[q], dict) and 'data' in data[q]:
            data2[q2] = data[q]['data']
        else:
            data2[q2] = data[q]
        if q in hasdope and q in hastemp:
            # temperature index first for consistency with other codes
            if q in hastype:
                for t in data2[q2]:
                    data2[q2][t] = np.swapaxes(data2[q2][t],0,1)
            else:
                data2[q2] = np.swapaxes(data2[q2],0,1)
        if q in hastype:
            if 'scattering_labels' not in data2:
                data2['scattering_labels'] = data[q].keys()
            # for consistency with the format in the mesh data
            data2[q2] = [data2[q2][l] for l in data2['scattering_labels']]
        if q2 in units:
            data2['meta']['units'][q2] = units[q2]

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], conversions[c])

    return data2

def amset_mesh(filename, quantities=['temperature', 'doping',
               'scattering_rates', 'scattering_labels'], spin='avg'):
    """Loads AMSET mesh data.

    Includes unit conversion and outputs units (see tp.settings).
    Swaps temperature and doping indices so temperature is first, for
    consistency with other codes.

    Arguments
    ---------

        filename : str
            filepath.

        quantites : dict, optional
            values to extract. Accepts AMSET keys, without spin
            channels, which are dealt with in the spin variable.
            Default: temperature, doping, scattering_rates,
            scattering_labels.

        spin : str, optional
            spin. Accepts up, down or avg. If avg and there is only one
            spin channel, selects that, else averages both. Default: avg.

    Returns
    -------

        dict
            extracted values.
    """

    import h5py

    # name conversions and abbreviations

    conversions = settings.amset_conversions()
    anames = settings.to_amset()
    tnames = settings.to_tp()
    units = settings.units()
    if isinstance(quantities, str): quantities = quantities.split()
    quantities = [anames[q] if q in anames else q for q in quantities]

    # list of quantities dependant on doping and temperature
    hasdope = ['fermi_levels', 'scattering_rates']
    hastemp = ['fermi_levels', 'scattering_rates']
    hastype = ['scattering_rates']
    hasspin = ['energies', 'vb_index', 'scattering_rates', 'velocities']

    # add dependant variables

    if 'doping' not in quantities:
        for q in quantities:
            if q in hasdope:
                quantities.append('doping')
                break
    if 'temperatures' not in quantities:
        for q in quantities:
            if q in hastemp:
                quantities.append('temperatures')
                break
    if 'scattering_labels' not in quantities:
        for q in quantities:
            if q in hastype:
                quantities.append('scattering_labels')
                break

    # load data

    data = h5py.File(filename, 'r')

    if spin in ['avg', 'average']:
        if 'energies_up' in data and 'energies_down' in data:
            spin = 'avg'
        elif 'energies_up' in data:
            spin = 'up'
        elif 'energies_down' in data:
            spin = 'down'
    data2 = {'meta': {'electronic_source': 'amset',
                      'spin':              spin,
                      'units':             {}}}
    for q in quantities:
        q2 = tnames[q] if q in tnames else q
        if q in hasspin:
            if spin == 'avg':
                data2[q2] = np.average([data['{}_up'.format(q)],
                                        data['{}_down'.format(q)]], axis=0)
            elif spin in ['up', 'down']:
                data2[q2] = data['{}_{}'.format(q, spin)]
            else:
                raise Exception('spin must be up or down or avg')
        else:
            data2[q2] = data[q]
        if q in hasdope and q in hastemp:
            # temperature in first index for consistency with other codes
            if q in hastype:
                data2[q2] = np.swapaxes(data2[q2],1,2)
            else:
                data2[q2] = np.swapaxes(data2[q2],0,1)
        if q2 in units:
            data2['meta']['units'][q2] = units[q2]

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], conversions[c])

    return data2

def boltztrap(filename, quantities=['temperature', 'doping', 'seebeck',
              'conductivity', 'electronic_thermal_conductivity'], doping='n'):
    """Loads BoltzTraP data from the tp boltztrap.hdf5 file.

    Includes unit conversion and outputs units (see tp.settings).

    Arguments
    ---------

        filename : str
            filepath.

        quantites : dict, optional
            values to extract. Accepts boltztrap.hdf5 keys.
            Default: temperature, doping, seebeck, conductivity,
            electronic_thermal_conductivity.

        doping : str, optional
            doping.  Default: n.

    Returns
    -------

        dict
            extracted values.
    """

    import h5py

    # name conversions and abbreviations

    assert doping in ['n', 'p'], 'doping must be n or p'

    conversions = settings.boltztrap_conversions()
    bnames = settings.to_boltztrap()
    tnames = settings.to_tp()
    units = settings.units()
    if isinstance(quantities, str): quantities = quantities.split()
    quantities = [bnames[q] if q in bnames else q for q in quantities]

    # list of quantities dependant on doping and temperature
    hasdope = ['average_eff_mass', 'conductivity', 'fermi_level', 'seebeck',
               'power_factor', 'electronic_thermal_conductivity']
    hastemp = ['average_eff_mass', 'conductivity', 'fermi_level', 'seebeck',
               'power_factor', 'electronic_thermal_conductivity']

    # add dependant variables

    if 'doping' not in quantities:
        for q in quantities:
            if q in hasdope:
                quantities.append('doping')
                break
    if 'temperature' not in quantities:
        for q in quantities:
            if q in hastemp:
                quantities.append('temperature')
                break

    # load data

    data = h5py.File(filename, 'r')

    data2 = {'meta': {'electronic_source': 'boltztrap',
                      'units':             {}}}
    for q in quantities:
        assert q in data, '{} unrecognised. Quantity must be in {} or {}.'.format(
                q, ', '.join(list(data)[:-1]), list(data)[-1])
        q2 = tnames[q] if q in tnames else q
        data2[q2] = data[q]
        if q in hasdope: data2[q2] = data2[q2][doping]
        if q2 in units:
            data2['meta']['units'][q2] = units[q2]

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], conversions[c])

    return data2

def phono3py(filename, quantities=['kappa', 'temperature'],
             write_lifetime=False, write_mfp=False, write_occupation=False):
    """Loads Phono3py data.

    Can also calculate lifetimes, mean free paths and occupations, which
    can be written to a file.
    Includes unit conversions and outputs units for all the data (see
    tp.settings). Also corrects mode_kappa for different phono3py
    versions.

    Arguments
    ---------

        filename : str
            filepath.

        quantities : list, optional
            values to extract. Accepts Phono3py keys, lifetime,
            mean_free_path and occupation. Default: kappa, temperature.
        write_lifetime : bool, optional
            write lifetimes to a new hdf5 file if in quantites.
            Default: False.
        write_mfp : bool, optional
            write mean free paths to a new hdf5 file if in quantities.
            Default: False.
        write_occupation : bool, optional
            write occupations to a new hdf5 file if in quantities.
            Default: False.

    Returns
    -------

        dict
            output data.
    """

    import h5py

    # name conversions and abbreviations

    conversions = settings.phono3py_conversions()
    pnames = settings.to_phono3py()
    tnames = settings.to_tp()
    units = settings.units()
    if isinstance(quantities, str): quantities = quantities.split()
    quantities = [pnames[q] if q in pnames else q for q in quantities]
    subs = {'dispersion': 'qpoint',
            'waterfall':  'frequency',
            'wideband':   ['frequency', 'gamma', 'qpoint']}
    hast = ['gamma', 'heat_capacity', 'kappa', 'lifetime',
            'mean_free_path','mode_kappa', 'occupation']
    for i in range(len(quantities)):
        if quantities[i] in subs:
            quantities[i] = subs[quantities[i]]

    quantities = list(np.ravel(quantities))

    # add dependant variables

    if 'temperature' not in quantities:
        for q in quantities:
            if q in hast:
                quantities.append('temperature')
                break

    # load and calculate data

    data = h5py.File(filename, 'r')

    data2 = {'meta': {'kappa_source': 'phono3py',
                      'units':        {}}}
    for q in quantities:
        assert q in data or q in ['lifetime', 'mean_free_path', 'occupation'], \
           '{} unrecognised. Quantity must be {}, lifetime, mean_free_path ' \
           'or occupation'.format(q, ', '.join(data))
        q2 = tnames[q] if q in tnames else q
        if q in data:
            data2[q2] = data[q][()]
        elif q in ['lifetime', 'mean_free_path']:
            data2['lifetime'] = np.reciprocal(np.multiply(2,data['gamma'][()]))
            data2['lifetime'] = np.where(np.isinf(data2['lifetime']), 0,
                                         data2['lifetime'])
            if q == 'mean_free_path':
                data2[q] = np.multiply(np.transpose([data2['lifetime'],] * 3,
                                                    (1, 2, 3, 0)),
                                       data['group_velocity'][()])
        elif q == 'occupation':
            from tp.calculate import be_occupation as occupation
            data2[q] = [occupation(data['frequency'][()], t)
                        for t in data['temperature'][()]]
        if q2 in units:
            data2['meta']['units'][q2] = units[q2]

    # write calculated data (loath to mess with original file)
    for write, q in zip([write_lifetime, write_mfp, write_occupation],
                        ['lifetime', 'mean_free_path', 'occupation']):
        if write and q in quantities:
            data3 = h5py.File('{}-{}'.format(q, filename), 'w')
            for q2 in [q, 'temperature', 'qpoint']:
                data3.create_dataset(q2, np.shape(data2[q2]), data=data2[q2])
            data3.close()

    # check mode_kappa and correct for certain phono3py versions
    if 'mode_kappa' in data2:
        try:
            k = round(data['kappa'][-1][0], 3)
            mk = round(data['mode_kappa'][-1][:,:,0].sum(axis=1).sum(), 3)
            if k != mk:
                raise Exception('The sum of mode_kappa does not equal kappa.\n '
                                'kappa={:.3f}; sum(mode_kappa)={:.3f}.'.format(
                                                                        k, mk))
        except Exception:
            mk2 = np.divide(data['mode_kappa'], np.prod(data['mesh'][:]))
            k = round(data['kappa'][-1][0], 3)
            mk = round(mk2[-1][:,:,0].sum(axis=1).sum(), 3)
            if k != mk:
                raise Exception('Mode kappa has been divided by the mesh, but '
                                'the sum of mode_kappa does not equal kappa.\n '
                                'kappa={:.3f}; sum(mode_kappa)={:.3f}.'.format(
                                                                        k, mk))
            else:
                data2['mode_kappa'] = np.divide(data2['mode_kappa'],
                                                np.prod(data['mesh'][()][:]))

    data.close()

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], conversions[c])

    return data2

def phonopy_dispersion(filename, xdata=None):
    """Loads phonopy dispersion, and can scale the x values.

    Scaling the x values is necessary to plot multiple dispersions on
    the same axes.

    Arguments
    ---------

        filename : str
            filepath.

        xdata : dict, optional
            data for the dispersion to scale this to. Should have the
            same path, must have the same number of labels. Not
            necessary if using tp.plot.phonons.add_multi. Default: None.

    Returns
    -------

        dict
            dispersion data.
    """

    import yaml

    # load data

    with open(filename, 'r') as f:
        data = yaml.safe_load(f)

    x = [d['distance'] for d in data['phonon']]
    qp = [q['q-position'] for q in data['phonon']]
    d2, ticks = get_path(data)
    eigs = [[b['frequency'] for b in p['band']] for p in data['phonon']]

    # scale data to other path

    if xdata is not None:
        d1 = xdata['tick_position']
        n = 0
        for i, d0 in enumerate(x):
            while n <= len(d2) and not (d0 >= d2[n] and d0 <= d2[n+1]):
                n += 1
            x[i] = d1[n] + ((d0 - d2[n]) * (d1[n+1] - d1[n]) / \
                                           (d2[n+1] - d2[n]))
    else:
        d1 = d2

    conversions = settings.phonopy_conversions()
    if 'frequency' in conversions:
        eigs *= conversions['frequency']

    units = tp.settings.units()
    data2 = {'x':             x,
             'qpoint':        qp,
             'frequency':     eigs,
             'tick_position': d1,
             'tick_label':    ticks,
             'meta':
                 {'phonon_dispersion_source': 'phonopy',
                     'units': {'frequency': units['frequency']}}}

    return data2

def phonopy_dos(filename, poscar='POSCAR', atoms=None):
    """Loads phonopy DoS and collates data per atom and as a total.

    By default reads atom names from a POSCAR, but can be overridden to
    allow for separation of environments.

    Arguments
    ---------

        filename : str
            path to phonopy projected_dos.dat or similar.
        poscar : str, optional
            path to POSCAR. Ignored if atoms specified. Default: POSCAR.
        atoms : str or array-like, optional
            atoms in POSCAR order. Atom names can be repeated, in which
            case their contributions are summed. Numbers can indicate
            repetition in the manner of a chemical formula, so the
            following are all acceptable and equivalent: "Ba 1 Sn 2 O 3",
            "Ba Sn Sn O O O", "Ba Sn O 3". Different environments can be
            distinguised with different atom names.
            Default: read from POSCAR.

    Returns
    -------

        dict
            frequency, DoS per atom and total.
    """

    # load data

    data = np.transpose(np.loadtxt(filename))
    units = tp.settings.units()
    data2 = {'frequency': data[0],
             'meta':      {'phonon_dos_source': 'phonopy',
                           'units': {'frequency': units['frequency']}}}

    conversions = settings.phonopy_conversions()
    if 'frequency' in conversions:
        data2['frequency'] *= conversions['frequency']

    if atoms is None:
        from pymatgen.io.vasp.inputs import Poscar
        poscar = Poscar.from_file(poscar, check_for_POTCAR=False,
                                  read_velocities=False).as_dict()
        atoms = [p['label'] for p in poscar['structure']['sites']]
    elif isinstance(atoms, str):
        atoms = atoms.split()

    # combine atoms contributions

    i = 0
    n = 1
    while i < len(atoms):
        try:
            atoms[i+1] = int(atoms[i+1])
            if atoms[i] in data2:
                data2[atoms[i]] += np.sum(data[n:n+atoms[i+1]], axis=0)
            else:
                data2[atoms[i]] = np.sum(data[n:n+atoms[i+1]], axis=0)
            n += atoms[i+1]
            i += 2
        except Exception:
            if atoms[i] in data2:
                data2[atoms[i]] += data[n]
            else:
                data2[atoms[i]] = data[n]
            n += 1
            i += 1

    data2['total'] = np.sum(data[1:], axis=0)

    return data2

def get_path(yamldata):
    """Extracts the path from a phonopy yaml.

    Works for the old and new phonopy and sumo yaml files.

    Arguments
    ---------

        yamldata : str
            raw phonopy dispersion data (i.e. from yaml.safe_load).

    Returns
    -------

        list
            x tick ordinates.
        list
            x tick labels.
    """

    tickindex = [1]
    tickindex.extend(np.cumsum(yamldata['segment_nqpoint']))
    tickpos = [yamldata['phonon'][i-1]['distance'] for i in tickindex]

    try: # old phonopy
        ticks = [i[0] for i in yamldata['labels']]
        ticks.append(yamldata['labels'][-1][1])
    except: # new phonopy/ sumo
        ticks = [yamldata['phonon'][i - 1]['label'] for i in tickindex]
    ticks = ['$\mathregular{\Gamma}$' if i == 'G' or 'gamma' in i.lower() else
             '$\mathregular{{{}}}$'.format(i.strip('$')) for i in ticks]

    return tickpos, ticks
