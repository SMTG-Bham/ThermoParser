"""Data loading tools.

Functions:
    amset
    phono3py
    phonopy_dispersion
    phonopy_dos

    get_path:
        gets high path from phonopy dispersion data.
"""

import numpy as np
import tp
from tp import settings
from . import aniso

def amset(filename, quantities=['temperatures', 'doping', 'seebeck',
          'conductivity', 'electronic_thermal_conductivity'], spin='up'):
    """Loads AMSET data.

    Includes unit conversion and outputs units (see tp.settings).

    Arguments:
        filename : str
            filepath.

        quantites : dict, optional
            values to extract. Accepts AMSET keys
            Default: temperatures, doping, seebeck, conductivity,
            electronic_thermal_conductivity.

        spin : str, optional
            spin.  Default: up.

    Returns:
        dict
            extracted values.
    """

    import json

    conversions = settings.conversions()
    anames = settings.to_amset()
    tnames = settings.to_tp()
    units = settings.units()
    quantities = [anames[q] if q in anames else q for q in quantities]

    # list of quantities requiring 'data' key, scattering type and spin
    hasdata = ['doping', 'temperatures', 'fermi_levels', 'conductivity',
               'seebeck', 'electronic_thermal_conductivity',# 'mobility',
               'kpoints', 'ir_kpoints', 'ir_to_full_kpoint_mapping']
    hasdope = ['fermi_levels', 'conductivity', 'seebeck',
               'electronic_thermal_conductivity', 'mobility', 'scattering_rates']
    hastemp = ['fermi_levels', 'conductivity', 'seebeck',
               'electronic_thermal_conductivity', 'mobility', 'energies',
               'velocities_product', 'scattering_rates']
    #hastype = ['mobility', 'scattering_rates', 'scattering_labels']
    hasspin = ['energies', 'vb_idx', 'velocities_product']#, 'scattering_rates']

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
        assert q in data, '{} unrecognised. Quantity must be in {}.'.format(
                           q, ', '.join(list(data)))
        q2 = tnames[q] if q in tnames else q
        data2[q2] = data[q]
        if q in hasspin: data2[q2] = data2[q2][spin]
        if q in hasdata: data2[q2] = data2[q2]['data']
        if q2 in units:
            data2['meta']['units'][q2] = units[q2]

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], conversions[c])

    return data2

def phono3py(filename, quantities=['kappa', 'temperature'], temperature=300,
             direction='norm'):
    """Loads Phono3py data.

    Can also calculate lifetimes, mean free paths. and occupations.
    Includes unit conversions and outputs units for all the data (see
    tp.settings). Also corrects mode_kappa for different phono3py
    versions.

    Arguments:
        filename : str
            filepath.

        quantities : list, optional
            values to extract. Accepts Phono3py keys, lifetime,
            mean_free_path and occupation. Default: kappa, temperature.

        temperature : float, optional
            temperature in K for mean_free_path. Default: 300.
        direction : str, optional
            direction for mean_free_path. Accepts a-c/ x-y,
            average/ avg or normal/ norm. Default: normal.

    Returns:
        dict
            output data.
    """

    import h5py

    conversions = settings.conversions()
    pnames = settings.to_phono3py()
    tnames = settings.to_tp()
    units = settings.units()
    if isinstance(quantities, str): quantities = quantities.split()
    quantities = [pnames[q] if q in pnames else q for q in quantities]
    subs = {'dispersion': 'qpoint',
            'waterfall':  'frequency'}
    hast = ['gamma', 'heat_capacity', 'kappa', 'lifetime', 'mode_kappa',
            'occupation']
    for i in range(len(quantities)):
        if quantities[i] in subs:
            quantities[i] = subs[quantities[i]]

    if 'temperature' not in quantities:
        for q in quantities:
            if q in hast:
                quantities.append('temperature')
                break

    # load data

    data = h5py.File(filename, 'r')
    ti = np.abs(np.subtract(data['temperature'], temperature)).argmin()

    data2 = {'meta': {'kappa_source': 'phono3py',
                      'units':        {}}}
    for q in quantities:
        assert q in data or q in ['lifetime', 'mean_free_path', 'occupation'], \
           '{} unrecognised. Quantity must be {}, lifetime, mean_free_path. ' \
           'or occupation'.format(q, ', '.join(data))
        q2 = tnames[q] if q in tnames else q
        if q in data:
            data2[q2] = data[q]
        elif q in ['lifetime', 'mean_free_path']:
            data2['lifetime'] = np.reciprocal(np.multiply(2,data['gamma']))
            lmax = np.amax(np.ma.masked_invalid(data2['lifetime']).compressed())
            for i in np.argwhere(np.isinf(data2['lifetime'])):
                data2['lifetime'][i[0],i[1]] = lmax
            if q == 'mean_free_path':
                data2[q] = np.multiply(data2['lifetime'][ti],
                                aniso.three(data['group_velocity'], direction))
        elif q == 'occupation':
            from tp.calculate import be_occupation as occupation
            data2[q] = [occupation(data['frequency'], t) for t in data['temperature']]

        if q2 in units:
            data2['meta']['units'][q2] = units[q2]

    # check mode_kappa and correct for certain phono3py versions

    if 'mode_kappa' in data2:
        try:
            k = round(data['kappa'][ti][0], 3)
            mk = round(data['mode_kappa'][ti][:,:,0].sum(axis=1).sum(), 3)
            if k != mk:
                raise Exception('The sum of mode_kappa does not equal kappa.\n '
                                'kappa={:.3f}; sum(mode_kappa)={:.3f}.'.format(
                                                                        k, mk))
        except Exception:
            mk2 = np.divide(data['mode_kappa'], np.prod(data['mesh'][:]))
            k = round(data['kappa'][ti][0], 3)
            mk = round(mk2[ti][:,:,0].sum(axis=1).sum(), 3)
            if k != mk:
                raise Exception('Mode kappa has been divided by the mesh, but '
                                'the sum of mode_kappa does not equal kappa.\n '
                                'kappa={:.3f}; sum(mode_kappa)={:.3f}.'.format(
                                                                        k, mk))
            else:
                data2['mode_kappa'] = np.divide(data2['mode_kappa'],
                                                np.prod(data['mesh'][:]))

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], conversions[c])

    return data2

def phonopy_dispersion(filename, xdata=None):
    """Loads phonopy dispersion, and can scale the x values.

    Scaling the x values is necessary to plot multiple dispersions on
    the same axes.

    Arguments:
        filename : str
            filepath.

        xdata : dict, optional
            data for the dispersion to scale this to. Should have the
            same path, must have the same number of labels. Default: None.

    Returns:
        dict
            dispersion data.
    """

    import yaml

    with open(filename, 'r') as f:
        data = yaml.safe_load(f)

    x = [d['distance'] for d in data['phonon']]
    qp = [q['q-position'] for q in data['phonon']]
    d2, ticks = get_path(data)
    eigs = [[b['frequency'] for b in p['band']] for p in data['phonon']]

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

    units = tp.settings.units()
    data2 = {'x':             x,
             'qpoint':        qp,
             'eigenvalue':    eigs,
             'tick_position': d1,
             'tick_label':    ticks,
             'meta':
                 {'phonon_dispersion_source': 'phonopy',
                  'units':
                      {'frequency': units['frequency']}}}

    return data2

def phonopy_dos(filename, atoms):
    """Loads phonopy DoS and collates data per atom and as a total.

    Atoms are user-defined, allowing separation of environments.

    Arguments:
        filename : str
            filepath.
        atoms : str or array-like
            atoms and numbers in the format "atom number" in POSCAR
            order, e.g. for BaSnO3: "Ba 1 Sn 2 O 3"

    Returns:
        dict
            DoS per atom and total.
    """

    data = np.transpose(np.loadtxt(filename))
    data2 = {'x':    data[0],
             'meta': {'phonon_dos_source': 'phonopy'}}

    if isinstance(atoms, str): atoms = atoms.split()

    i = 0
    n = 1
    while i < len(atoms):
        atoms[i+1] = int(atoms[i+1])
        data2[atoms[i]] = np.sum(data[n:n+atoms[i+1]], axis=0)
        n += atoms[i+1]
        i += 2
    data2['total'] = np.sum(data, axis=0)

    return data2

def get_path(yamldata):
    """Extracts the path from a phonopy yaml.

    Works for the old and new phonopy and sumo yaml files.

    Arguments:
        yamldata : str
            raw phonopy dispersion data (i.e. from yaml.safe_load).

    Returns:
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
    ticks = ['$\mathregular{\Gamma}$' if i == 'G' or i == '\Gamma' else
             '$\mathregular{{{}}}$'.format(i) for i in ticks]

    return tickpos, ticks
