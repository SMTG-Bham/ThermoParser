"""Data loading tools.

Loads data from codes into a dictionary, with units and array structures
standardised. Also adds a ``meta`` subdictionary, which contains units,
array dimensions and the data source.

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

def amset(filename, quantities=['seebeck', 'conductivity',
          'electronic_thermal_conductivity'], doping='n'):
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

        quantites : str or list, optional
            values to extract. Accepts AMSET keys, without spin
            channels, which are dealt with in the spin variable. Loads
            dependant properties. Default: seebeck, conductivity,
            electronic_thermal_conductivity.

        doping : str, optional
            doing type (n or p). If there is more than one, defaults to
            n, else this is ignored.

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
    dimensions = settings.dimensions()
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

    if 'doping' in quantities:
        d = np.array(data['doping'])
        if (d < 0).all():
            doping = 'n'
            di = np.where(d < 0)[0]
        elif (d > 0).all():
            doping = 'p'
            di = np.where(d > 0)[0]
        elif doping == 'n':
            di = np.where(d < 0)[0]
        elif doping == 'p':
            di = np.where(d > 0)[0]

    data2 = {'meta': {'doping_type':       doping,
                      'electronic_source': 'amset',
                      'units':             {},
                      'dimensions':        {}}}
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
        if q in hasdope:
            # temperature index first for consistency with other codes
            if q in hastype:
                for t in data2[q2]:
                    data2[q2][t] = np.array(data2[q2][t])[di]
                    if q in hastemp:
                        data2[q2][t] = np.swapaxes(data2[q2][t],0,1)
            else:
                data2[q2] = np.array(data2[q2])[di]
                if q in hastemp:
                    data2[q2] = np.swapaxes(data2[q2],0,1)
        if q in hastype:
            if 'scattering_labels' not in data2:
                data2['scattering_labels'] = data[q].keys()
            # for consistency with the format in the mesh data
            data2[q2] = [data2[q2][l] for l in data2['scattering_labels']]
        if q2 in units:
            data2['meta']['units'][q2] = units[q2]
        if q2 in dimensions:
            data2['meta']['dimensions'][q2] = dimensions[q2]

    if 'doping' in data2:
        data2['doping'] = np.array(data2['doping'])[di]

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], conversions[c])

    return data2

def amset_mesh(filename, quantities='scattering_rates', doping='n',
               spin='avg'):
    """Loads AMSET mesh data.

    Can also weight rates. Includes unit conversion and outputs units
    (see tp.settings). Swaps temperature and doping indices so
    temperature is first, for consistency with other codes.

    Arguments
    ---------

        filename : str
            filepath.

        quantites : str or list, optional
            values to extract. Accepts AMSET keys, without spin
            channels, which are dealt with in the spin variable.
            Also accepts ibz_weights, the weights of the irreducible
            k-points, fd_weights, the weights of the energies wrt the
            derivative of the Fermi-Dirac distribution, and
            weighted_rates, scattering_rates weighted by fd_weights
            and averaged over kpoints. Loads dependant properties.
            Default: scattering_rates.

        doping : str, optional
            doing type (n or p). If there is more than one, defaults to
            n, else this is ignored.

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
    dimensions = settings.dimensions()
    if isinstance(quantities, str): quantities = quantities.split()
    quantities = [anames[q] if q in anames else q for q in quantities]

    # list of abbriviations and dependant quantites
    subs = {'weights': ['ibz_weights', 'fd_weights']}
    hasdope = ['fd_weights', 'fermi_levels', 'normalised_weights',
               'scattering_rates', 'weighted_rates']
    hastemp = ['fd_weights', 'energies', 'fermi_levels', 'normalised_weights',
               'scattering_rates', 'weighted_rates']
    hastype = ['scattering_rates', 'weighted_rates']
    hasspin = ['energies', 'vb_index', 'scattering_rates', 'velocities']

    for i in range(len(quantities)):
        if quantities[i] in subs:
            quantities[i] = subs[quantities[i]]

    quantities = list(np.ravel(quantities))

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

    def resolve_spin(data, q, spin):
        """Resolves spin in AMSET.

        Either returns the up or down value, or averages them.

        Arguments
        ---------

        data : dict
            data dictionary from AMSET.
        q : str
            quantity to resolve
        spin : str
            spin to use. Must be up, down or average.

        Returns
        -------

        array
            resolved data
        """

        if spin in ['avg', 'average']:
            resolved = np.average([data['{}_up'.format(q)][()],
                                   data['{}_down'.format(q)][()]], axis=0)
        elif spin in ['up', 'down']:
            resolved = data['{}_{}'.format(q, spin)][()]
        else:
            raise Exception('spin must be up or down or average')

        return resolved

    with h5py.File(filename, 'r') as f:
        if 'doping' in quantities:
            d = np.array(f['doping'][()])
            if (d < 0).all():
                doping = 'n'
                di = np.where(d < 0)[0]
            elif (d > 0).all():
                doping = 'p'
                di = np.where(d > 0)[0]
            elif doping == 'n':
                di = np.where(d < 0)[0]
            elif doping == 'p':
                di = np.where(d > 0)[0]

        if spin in ['avg', 'average']:
            if 'energies_up' in f and 'energies_down' in f:
                spin = 'avg'
            elif 'energies_up' in f:
                spin = 'up'
            elif 'energies_down' in f:
                spin = 'down'

        data = {'meta': {'doping_type':       doping,
                         'electronic_source': 'amset',
                         'spin':              spin,
                         'units':             {},
                         'dimensions':        {}}}
        for q in quantities:
            q2 = tnames[q] if q in tnames else q
            if q in hasspin:
                data[q2] = resolve_spin(f, q, spin)
            elif q in f:
                data[q2] = f[q][()]
            elif q not in ['ibz_weights', 'fd_weights', 'weighted_rates']:
                raise Exception('{} unrecognised. Quantity must be {}, '
                                'ibz_weights, fd_weights weighted_rates'
                                ''.format(q, ', '.join(f)))
            if q in ['ibz_weights', 'weighted_rates']:
                _, data['ibz_weights'] = np.unique(f['ir_to_full_kpoint_mapping'],
                                                   return_counts=True)
            if q in ['fd_weights', 'weighted_rates']:
                e = resolve_spin(f, 'energies', spin)
                data['fd_weights'] = -tp.calculate.dfdde(e, f['fermi_levels'],
                                                         f['temperatures'],
                                                         f['doping'],
                                                         amset_order=True)
            if q == 'weighted_rates':
                rates = resolve_spin(f, 'scattering_rates', spin)
                rates[rates > 1e20] = 1e15
                data['total_weights'] = data['fd_weights'] * data['ibz_weights']
                data['normalised_weights'] = data['total_weights'] / \
                                             np.sum(data['total_weights'],
                                                    axis=(2,3))[:,:,None,None]
                data[q2] = rates * data['normalised_weights']
                data[q2] = data[q2].sum(axis=(3,4))

        for q2 in data:
            q = anames[q2] if q2 in anames else q2
            if q in hasdope:
                # temperature in first index for consistency with other codes
                if q in hastype:
                    if q in hastemp:
                        data[q2] = np.swapaxes(data[q2],1,2)
                        data[q2] = np.array(data[q2])[:,:,di]
                    else:
                        data[q2] = np.array(data[q2])[:,di]
                else:
                    if q in hastemp:
                        data[q2] = np.swapaxes(data[q2],0,1)
                        data[q2] = np.array(data[q2])[:,di]
                    else:
                        data[q2] = np.array(data[q2])[di]
            if q2 in units:
                data['meta']['units'][q2] = units[q2]
            if q2 in dimensions:
                data['meta']['dimensions'][q2] = dimensions[q2]

    if 'doping' in data:
        data['doping'] = np.array(data['doping'])[di]

    if 'scattering_labels' in data:
        data['scattering_labels'] = \
                         [l.decode('UTF-8') for l in data['scattering_labels']]

    for c in conversions:
        if c in data:
            data[c] = np.multiply(data[c], conversions[c])

    return data

def boltztrap(filename, quantities=['temperature', 'doping', 'seebeck',
              'conductivity', 'electronic_thermal_conductivity'], doping='n'):
    """Loads BoltzTraP data from the tp boltztrap.hdf5 file.

    Includes unit conversion and outputs units (see tp.settings).

    Arguments
    ---------

        filename : str
            filepath.

        quantites : str or list, optional
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
    dimensions = settings.dimensions()
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

    with h5py.File(filename, 'r') as f:
        data = {'meta': {'doping_type':       doping,
                         'electronic_source': 'boltztrap',
                         'units':             {},
                         'dimensions':        {}}}
        for q in quantities:
            assert q in f, '{} unrecognised. Quantity must be in {} or {}.'.format(
                            q, ', '.join(list(f)[:-1]), list(f)[-1])
            q2 = tnames[q] if q in tnames else q
            if q in hasdope:
                data[q2] = f[q][doping][()]
            else:
                data[q2] = f[q][()]
            if q2 in units:
                data['meta']['units'][q2] = units[q2]
            if q2 in dimensions:
                data['meta']['dimensions'][q2] = dimensions[q2]

    for c in conversions:
        if c in data:
            data[c] = np.multiply(data[c], conversions[c])

    return data

def phono3py(filename, quantities=['kappa', 'temperature']):
    """Loads Phono3py data.

    Can also calculate lifetimes, mean free paths and occupations.
    Includes unit conversions and outputs units for all the data (see
    tp.settings). Also corrects mode_kappa for different phono3py
    versions.

    Arguments
    ---------

        filename : str
            filepath.

        quantities : str or list, optional
            values to extract. Accepts Phono3py keys, lifetime,
            mean_free_path and occupation. Default: kappa, temperature.

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
    dimensions = settings.dimensions()
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
            if q in hast or (q in tnames and tnames[q] in hast):
                quantities.append('temperature')
                break

    # load and calculate data

    with h5py.File(filename, 'r') as f:
        data = {'meta': {'kappa_source': 'phono3py',
                         'units':        {},
                         'dimensions':   {}}}
        for q in quantities:
            assert q in f or q in ['lifetime', 'mean_free_path', 'occupation'], \
               '{} unrecognised. Quantity must be {}, lifetime, mean_free_path ' \
               'or occupation'.format(q, ', '.join(f))
            q2 = tnames[q] if q in tnames else q
            if q in f:
                data[q2] = f[q][()]
            elif q in ['lifetime', 'mean_free_path']:
                data['lifetime'] = np.reciprocal(np.multiply(2 * 2 * np.pi,
                                                             f['gamma'][()]))
                data['lifetime'] = np.where(np.isinf(data['lifetime']), 0,
                                            data['lifetime'])
                if q == 'mean_free_path':
                    data[q] = np.multiply(np.transpose([data['lifetime'],] * 3,
                                                        (1, 2, 3, 0)),
                                          f['group_velocity'][()])
            elif q == 'occupation':
                from tp.calculate import be_occupation as occupation
                data[q] = [occupation(f['frequency'][()], t)
                            for t in f['temperature'][()]]
            if q2 in units:
                data['meta']['units'][q2] = units[q2]
            if q2 in dimensions:
                data['meta']['dimensions'][q2] = dimensions[q2]

        # check mode_kappa and correct for certain phono3py versions
        if 'mode_kappa' in data:
            try:
                k = round(f['kappa'][()][-1][0], 3)
                mk = round(f['mode_kappa'][()][-1][:,:,0].sum(axis=1).sum(), 3)
                if k != mk:
                    raise Exception('The sum of mode_kappa does not equal kappa.\n '
                                    'kappa={:.3f}; sum(mode_kappa)={:.3f}.'.format(
                                                                            k, mk))
            except Exception:
                mk2 = np.divide(f['mode_kappa'][()], np.prod(f['mesh'][()][:]))
                k = round(f['kappa'][()][-1][0], 3)
                mk = round(mk2[-1][:,:,0].sum(axis=1).sum(), 3)
                if k != mk:
                    raise Exception('Mode kappa has been divided by the mesh, but '
                                    'the sum of mode_kappa does not equal kappa.\n '
                                    'kappa={:.3f}; sum(mode_kappa)={:.3f}.'.format(
                                                                            k, mk))
                else:
                    data['mode_kappa'] = np.divide(data['mode_kappa'],
                                                   np.prod(f['mesh'][()][:]))

    for c in conversions:
        if c in data:
            data[c] = np.multiply(data[c], conversions[c])

    return data

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
    dimensions = settings.dimensions()
    data2 = {'x':             x,
             'qpoint':        qp,
             'frequency':     eigs,
             'tick_position': d1,
             'tick_label':    ticks,
             'meta':
                 {'phonon_dispersion_source': 'phonopy',
                     'units':      {'frequency': units['frequency']},
                     'dimensions': {'frequency': dimensions['frequency']}}}

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
    dimensions = settings.dimensions()
    data2 = {'frequency': data[0],
             'meta':      {'phonon_dos_source': 'phonopy',
                           'units':      {'frequency': units['frequency']},
                           'dimensions': {'frequency': dimensions['frequency']}}}

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

        yamldata : dict
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
