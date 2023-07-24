"""Data loading tools.

Loads data from codes into a dictionary, with units and array structures
standardised. Also adds a ``meta`` subdictionary, which contains units,
array dimensions and the data source.

Functions
---------

    amset:
        from the amset transport json.
    amset_mesh:
        from the amset mesh h5.
    boltztrap:
        from the tp boltztrap hdf5.
    phono3py:
        from the phono3py kappa hdf5.
    phonopy_dispersion:
        from the phonopy or sumo band.yaml.
    phonopy_dos:
        from the phonopy projected_dos.dat.

    get_path:
        gets high path from phonopy dispersion data.
"""

import numpy as np
import tp
from tp import settings

def amset(filename, quantities='all', doping='n'):
    """Loads AMSET transport data from json.

    Includes unit conversion and outputs units (see tp.settings).
    Swaps temperature and doping indices so temperature is first, for
    consistency with other codes, and mobility is formatted as an array,
    like scattering_rates is in the mesh.h5 file. Maintains basic
    compatibility with amset 0.1.

    Arguments
    ---------

        filename : str
            amset transport json filepath.

        quantites : str or list, optional
            values to extract. Accepts AMSET keys and power_factor.
            All loads the whole file, power factor can be added if it
            is not in the file. Loads dependent properties.
            Default: all.
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

    aconversions = settings.amset_conversions()
    conversions = settings.conversions()
    anames = settings.to_amset()
    tnames = settings.to_tp()
    units = settings.units()
    dimensions = settings.dimensions()
    if isinstance(quantities, str): quantities = quantities.split()
    quantities = [anames[q] if q in anames else q for q in quantities]

    # quantities derived from those in the file
    derived = {'power_factor': ['conductivity', 'seebeck']}

    # add dependent variables
    for d in ['doping', 'temperatures']:
        if d not in quantities:
            for q in quantities:
                if q in tnames:
                    q = tnames[q]
                if q in dimensions and \
                   (d in dimensions[q] or (d in tnames and tnames[d] in dimensions[q])):
                    quantities.append(d)
                    break
    for q in derived:
        if q in quantities:
            for q2 in derived[q]:
                if q2 not in quantities:
                    quantities.append(q2)

    # load data

    with open(filename) as f:
        data = json.load(f)
    if 'all' in quantities:
        if 'power_factor' in quantities:
            quantities = [*data.keys(), 'power_factor']
        else:
            quantities = list(data.keys())

    if 'doping' in quantities:
        if 'data' in data['doping']:
            d = np.array(data['doping']['data'])
        else:
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
        if q in derived:
            continue
        qs = [*list(data), *list(derived)]
        assert q in data, \
               '{} unrecognised. Quantity must be in {} or {}.'.format(q,
               ', '.join(qs[:-1]), qs[-1])
        q2 = tnames[q] if q in tnames else q
        # compatibility with previous version
        if isinstance(data[q], dict) and 'data' in data[q]:
            data2[q2] = data[q]['data']
        else:
            data2[q2] = data[q]
        if q2 in dimensions and 'doping' in dimensions[q2]:
            # temperature index first for consistency with other codes
            # With the latest version of resolve, this is unneccessary 
            # Should it be removed?
            if 'stype' in dimensions[q2]:
                for t in data2[q2]:
                    data2[q2][t] = np.array(data2[q2][t])[di]
                    if 'temperature' in dimensions[q2]:
                        data2[q2][t] = np.swapaxes(data2[q2][t],0,1)
            else:
                data2[q2] = np.array(data2[q2])[di]
                if 'temperature' in dimensions[q2]:
                    data2[q2] = np.swapaxes(data2[q2],0,1)
        if q2 in dimensions and 'stype' in dimensions[q2]:
            if q2 == 'mobility':
                data[q2]['Total'] = data[q2].pop('overall')
            if 'stype' not in data2:
                data2['stype'] = list(data[q].keys())
            # for consistency with the format in the mesh data
            data2[q2] = [data2[q2][l] for l in data2['stype']]
        if q2 in units:
            data2['meta']['units'][q2] = units[q2]
        if q2 in dimensions:
            data2['meta']['dimensions'][q2] = dimensions[q2]

    if 'doping' in data2:
        data2['doping'] = np.array(data2['doping'])[di]

    for c in aconversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], float(aconversions[c]))

    if 'power_factor' in quantities:
        data2 = tp.calculate.power_factor_fromdict(data2)

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], float(conversions[c]))

    return data2

def amset_mesh(filename, quantities='all', doping='n', spin='avg'):
    """Loads AMSET mesh data from h5.

    Can also weight rates. Includes unit conversion and outputs units
    (see tp.settings). Swaps temperature and doping indices so
    temperature is first, for consistency with other codes.

    Arguments
    ---------

        filename : str
            amset mesh h5 filepath.

        quantites : str or list, optional
            values to extract. Accepts AMSET keys, without spin
            channels, which are dealt with in the spin variable.
            Also accepts ibz_weights, the weights of the irreducible
            k-points, fd_weights, the weights of the energies wrt the
            derivative of the Fermi-Dirac distribution, weighted_rates,
            scattering_rates weighted by fd_weights and averaged over
            kpoints and occupation, the Fermi-Dirac occupation.
            "all" loads all quantities in the file, which does not
            include ibz_weights, fd_weights, weighted_rates or
            occupation. Loads dependent properties. Default: all.

        doping : str, optional
            doing type (n or p). If there is more than one, defaults to
            n, else this is ignored.

        spin : str, optional
            spin. Accepts up, down or avg. If avg and there is only one
            spin channel, selects that, else averages both. If avg,
            assumes non-spin separated bands and so multiplies occupation
            by 2. Default: avg.

    Returns
    -------

        dict
            extracted values.
    """

    import h5py

    # name conversions and abbreviations

    aconversions = settings.amset_conversions()
    conversions = settings.conversions()
    anames = settings.to_amset()
    tnames = settings.to_tp()
    units = settings.units()
    dimensions = settings.dimensions()
    dimensions['occupation'] = ['temperature', 'doping', 'kpoint', 'band']
    if isinstance(quantities, str): quantities = quantities.split()
    quantities = [anames[q] if q in anames else q for q in quantities]

    # list of abbriviations and dependent quantites
    subs = {'weights':        ['ibz_weights', 'fd_weights'],
            'weighted_rates': ['scattering_rates', 'energies', 'fermi_levels', 'weighted_rates'],
            'mean_free_path': ['scattering_rates', 'velocities', 'mean_free_path'],
            'weighted_mfp':   ['scattering_rates', 'velocities', 'mean_free_path', 'energies', 'fermi_levels', 'weighted_mfp'],
            'occupation':     ['energies', 'fermi_levels', 'occupation']}
    hasspin = ['energies', 'vb_index', 'scattering_rates', 'velocities']

    l = len(quantities)
    i = 0
    while i < l:
        if quantities[i] in subs:
            quantities.extend(subs[quantities[i]])
            del quantities[i]
            l -= 1
        else:
            i += 1

    # add dependent variables

    for d in ['doping', 'ir_kpoints', 'temperatures', 'scattering_labels']:
        if d not in quantities:
            for q in quantities:
                if q in tnames:
                    q = tnames[q]
                if q in dimensions and \
                   (d in dimensions[q] or (d in tnames and tnames[d] in dimensions[q])):
                    quantities.append(d)
                    break

    if 'temperatures' in quantities:
        quantities.insert(0, quantities.pop(quantities.index('temperatures')))

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
        if 'all' in quantities:
            q2 = quantities
            quantities = list(f.keys())
            for q in ['ibz_weights', 'fd_weights', 'weighted_rates', 'occupation']:
                if q in q2:
                    quantities.append(q)
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

        # spin2 so if avg selected doubles occupation
        if spin in ['avg', 'average']:
            if 'energies_up' in f and 'energies_down' in f:
                spin2 = 'avg'
            elif 'energies_up' in f:
                spin2 = 'up'
            elif 'energies_down' in f:
                spin2 = 'down'

        data = {'meta': {'doping_type':       doping,
                         'electronic_source': 'amset',
                         'spin':              spin2,
                         'units':             {},
                         'dimensions':        {}}}
        for q in quantities:
            q2 = tnames[q] if q in tnames else q
            if q in hasspin:
                data[q2] = resolve_spin(f, q, spin2)
            elif q in f:
                data[q2] = f[q][()]
            elif q not in ['ibz_weights', 'fd_weights', 'weighted_rates',
                           'mean_free_path', 'weighted_mfp', 'occupation']:
                raise Exception('{} unrecognised. Quantity must be {}, '
                                'ibz_weights, fd_weights, weighted_rates, '
                                'mean_free_path or weighted_mfp or occupation.'
                                ''.format(q, ', '.join(f)))
            if q == 'scattering_rates':
                data[q2] = [*list(data[q2]), list(np.sum(data[q2], axis=0))]
            if q in ['ibz_weights', 'weighted_rates', 'weighted_mfp']:
                _, data['ibz_weights'] = np.unique(f['ir_to_full_kpoint_mapping'],
                                                   return_counts=True)
            if q in ['fd_weights', 'weighted_rates', 'weighted_mfp']:
                e = resolve_spin(f, 'energies', spin2)
                data['fd_weights'] = -tp.calculate.dfdde(e, f['fermi_levels'],
                                                         f['temperatures'],
                                                         f['doping'],
                                                         amset_order=True)
            if q in ['mean_free_path', 'weighted_mfp']:
                data['velocities'] = np.array(data['velocities'])
                data['scattering_rates'] = np.array(data['scattering_rates'])
                data['mean_free_path'] = data['velocities'][None, None, None, :, :, :] / \
                                         data['scattering_rates'][:, :, :, :, :, None]
            if q in ['weighted_rates', 'weighted_mfp']:
                data['total_weights'] = data['fd_weights'] * data['ibz_weights']
                data['normalised_weights'] = data['total_weights'] / \
                                             np.sum(data['total_weights'],
                                                    axis=(2,3))[:,:,None,None]
            if q =='weighted_rates':
                from copy import deepcopy
                rates = np.array(deepcopy(data['scattering_rates']))
                rates[rates > 1e20] = 1e15
                data[q2] = rates * data['normalised_weights']
                data[q2] = data[q2].sum(axis=(3,4))
            if q == 'weighted_mfp':
                data[q2] = data['mean_free_path'] * \
                           data['normalised_weights'][None, :, :, :, :, None]
                data[q2] = data[q2].sum(axis=(3,4))
            if q == 'occupation':
                data['occupation'] = tp.calculate.fd_occupation(
                                     data['energy'], data['temperature'],
                                     data['fermi_level'])
                if spin in ['avg', 'average']:
                    data['occupation'] *= 2

        for q2 in data:
            q = anames[q2] if q2 in anames else q2
            if q2 in dimensions and 'doping' in dimensions[q2]:
                # temperature in first index for consistency with other codes
                if 'stype' in dimensions[q2]:
                    if 'temperature' in dimensions[q2]:
                        data[q2] = np.swapaxes(data[q2],1,2)
                        data[q2] = np.array(data[q2])[:,:,di]
                    else:
                        data[q2] = np.array(data[q2])[:,di]
                else:
                    if 'temperature' in dimensions[q2]:
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

    if 'stype' in data:
        data['stype'] = [l.decode('UTF-8') for l in data['stype']]
        data['stype'].append('Total')

    for c in aconversions:
        if c in data:
            data[c] = np.multiply(data[c], float(aconversions[c]))

    for c in conversions:
        if c in data:
            data[c] = np.multiply(data[c], float(conversions[c]))

    return data

def boltztrap(filename, quantities='all', doping='n'):
    """Loads BoltzTraP data from the tp boltztrap.hdf5 file.

    Includes unit conversion and outputs units (see tp.settings).

    Arguments
    ---------

        filename : str
            boltztrap.hdf5 filepath.

        quantites : str or list, optional
            values to extract. Accepts boltztrap.hdf5 keys or all.
            Default: all.

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

    bconversions = settings.boltztrap_conversions()
    conversions = settings.conversions()
    bnames = settings.to_boltztrap()
    tnames = settings.to_tp()
    units = settings.units()
    dimensions = settings.boltztrap_dimensions()
    if isinstance(quantities, str): quantities = quantities.split()
    if 'all' not in quantities:
        quantities = [bnames[q] if q in bnames else q for q in quantities]

        # add dependent variables

        for d in ['doping', 'temperature']:
            if d not in quantities:
                for q in quantities:
                    if q in tnames:
                        q = tnames[q]
                    if q in dimensions and \
                       (d in dimensions[q] or (d in tnames and tnames[d] in dimensions[q])):
                        quantities.append(d)
                        break

    # load data

    with h5py.File(filename, 'r') as f:
        if 'all' in quantities:
            quantities = list(f.keys())
            quantities.remove('meta')
        data = {'meta': {'doping_type':       doping,
                         'electronic_source': 'boltztrap',
                         'units':             {},
                         'dimensions':        {}}}
        for q in quantities:
            assert q in f, '{} unrecognised. Quantity must be in {} or {}.'.format(
                            q, ', '.join(list(f)[:-1]), list(f)[-1])
            q2 = tnames[q] if q in tnames else q
            if q2 in dimensions and 'dtype' in dimensions[q2]:
                data[q2] = f[q][doping][()]
                dimensions[q2].remove('dtype')
            else:
                data[q2] = f[q][()]
            if q2 in units:
                data['meta']['units'][q2] = units[q2]
            if q2 in dimensions:
                data['meta']['dimensions'][q2] = dimensions[q2]

    for c in bconversions:
        if c in data:
            data[c] = np.multiply(data[c], float(bconversions[c]))

    for c in conversions:
        if c in data:
            data[c] = np.multiply(data[c], float(conversions[c]))

    return data

def phono3py(filename, quantities='all'):
    """Loads Phono3py data from kappa hdf5.

    Can also calculate lifetimes, mean free paths and occupations.
    Includes unit conversions and outputs units and index order for all
    the data (see tp.settings). Also corrects mode_kappa for different
    phono3py versions. Also converts the default 6x1 direction matrices
    into 3x3 ones for compatability with other codes.

    Arguments
    ---------

        filename : str
            phono3py kappa hdf5 filepath.

        quantities : str or list, optional
            values to extract. Accepts Phono3py keys, lifetime,
            mean_free_path and occupation. 'all' loads all Phono3py
            keys, but lifetime etc. must be loaded separately.
            Default: all.

    Returns
    -------

        dict
            output data.
    """

    import h5py

    # name conversions and abbreviations

    pconversions = settings.phono3py_conversions()
    conversions = settings.conversions()
    pnames = settings.to_phono3py()
    tnames = settings.to_tp()
    units = settings.units()
    dimensions = settings.dimensions()
    if isinstance(quantities, str): quantities = quantities.split()
    quantities = [pnames[q] if q in pnames else q for q in quantities]
    subs = {'dispersion': ['qpoint'],
            'waterfall':  ['frequency'],
            'wideband':   ['frequency', 'gamma', 'qpoint']}

    # quantities derived from those in the file
    derived = {'lifetime': ['gamma'],
               'mean_free_path': ['gamma', 'group_velocity'],
               'occupation': ['frequency']}

    for i in range(len(quantities)):
        while True:
            if quantities[i] in subs:
                quantities.extend(subs[quantities[i]])
                del quantities[i]
            else:
                break
        
    for q in derived:
        if q in quantities:
            for q2 in derived[q]:
                if q2 not in quantities:
                    quantities.append(q2)

    quantities = list(np.ravel(quantities))

    # add dependent variables

    for d in ['temperature', 'qpoint']:
        if d not in quantities:
            for q in quantities:
                if q in tnames:
                    q = tnames[q]
                if q in dimensions and d in dimensions[q]:
                    quantities.append(d)
                    break

    # load and calculate data

    with h5py.File(filename, 'r') as f:
        if 'all' in quantities:
            q2 = quantities
            quantities = list(f.keys())
            for q in derived:
                if q in q2:
                    quantities.append(q)
        data = {'meta': {'kappa_source': 'phono3py',
                         'units':        {},
                         'dimensions':   {}}}
        for q in quantities:
            q2 = tnames[q] if q in tnames else q
            if q2 in units:
                data['meta']['units'][q2] = units[q2]
            if q2 in dimensions:
                data['meta']['dimensions'][q2] = dimensions[q2]
            if q in derived:
                continue
            qs = [*list(f), *list(derived)]
            assert q in f, \
                   '{} unrecognised. Quantity must be in {} or {}.'.format(q,
                   ', '.join(qs[:-1]), qs[-1])
            data[q2] = f[q][()]
            if q2 in dimensions:
                for i, d in enumerate(dimensions[q2]):
                    if d == 6:
                        data[q2] = np.moveaxis(data[q2], i, 0)
                        data[q2] = [[data[q2][0], data[q2][5], data[q2][4]],
                                    [data[q2][5], data[q2][1], data[q2][3]],
                                    [data[q2][4], data[q2][3], data[q2][2]]]
                        data[q2] = np.moveaxis(data[q2], 0, i+1)
                        data[q2] = np.moveaxis(data[q2], 0, i+1)
                        dimensions[q2][i] = 3
                        dimensions[q2].insert(i, 3)

        # check mode_kappa and correct for certain phono3py versions
        if 'mode_kappa' in data:
            try:
                k = round(f['kappa'][()][-1][0], 3)
                mk = round(f['mode_kappa'][()][-1][:,:,0].sum(axis=1).sum(), 3)
                if k != mk:
                    raise ValueError('The sum of mode_kappa does not equal kappa.\n '
                                     'kappa={:.3f}; sum(mode_kappa)={:.3f}.'.format(
                                                                             k, mk))
            except ValueError:
                mk2 = np.divide(f['mode_kappa'][()], np.prod(f['mesh'][()][:]))
                k = round(f['kappa'][()][-1][0], 3)
                mk = round(mk2[-1][:,:,0].sum(axis=1).sum(), 3)
                if k != mk:
                    raise ValueError('Mode kappa has been divided by the mesh, but '
                                     'the sum of mode_kappa does not equal kappa.\n '
                                     'kappa={:.3f}; sum(mode_kappa)={:.3f}.'.format(
                                                                             k, mk))
                else:
                    data['mode_kappa'] = np.divide(data['mode_kappa'],
                                                   np.prod(f['mesh'][()][:]))

    for c in pconversions:
        if c in data:
            data[c] = np.multiply(data[c], float(pconversions[c]))

    if 'lifetime' in quantities or 'mean_free_path' in quantities:
        data['lifetime'] = tp.calculate.lifetime(data['gamma'], use_tprc=False)
        if 'mean_free_path' in quantities:
            data['mean_free_path'] = tp.calculate.mfp(data['gamma'],
                                                      data['group_velocity'],
                                                      use_tprc=False)
    if 'occupation' in quantities:
        data['occupation'] = tp.calculate.be_occupation(data['frequency'],
                                                        data['temperature'],
                                                        use_tprc=False)

    for c in conversions:
        if c in data:
            data[c] = np.multiply(data[c], float(conversions[c]))

    return data

def phonopy_dispersion(filename, xdata=None):
    """Loads phonopy dispersion, and can scale the x values.

    Scaling the x values is necessary to plot multiple dispersions on
    the same axes.

    Arguments
    ---------

        filename : str
            phonopy or sumo band.yaml filepath.

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

    pconversions = settings.phonopy_conversions()
    conversions = settings.conversions()

    # load data

    with open(filename, 'r') as f:
        data = yaml.safe_load(f)

    x = [d['distance'] for d in data['phonon']]
    qp = [q['q-position'] for q in data['phonon']]
    tickpos, ticks = get_path(data)
    f = [[b['frequency'] for b in p['band']] for p in data['phonon']]

    if xdata is not None:
        # scale data to other path
        x = scale_to_path(x, tickpos, xdata['tick_position'])
        tickpos = xdata['tick_position']

    units = tp.settings.units()
    dimensions = settings.dimensions()
    data2 = {'x':             x,
             'qpoint':        qp,
             'frequency':     f,
             'tick_position': tickpos,
             'tick_label':    ticks,
             'meta':
                 {'phonon_dispersion_source': 'phonopy',
                  'units':      {'frequency': units['frequency']},
                  'dimensions': {'frequency': dimensions['frequency']}}}

    for c in pconversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], float(pconversions[c]))

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], float(conversions[c]))

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
    pconversions = settings.phonopy_conversions()
    conversions = settings.conversions()
    units = tp.settings.units()
    dimensions = settings.dimensions()
    data2 = {'frequency': data[0],
             'meta':      {'phonon_dos_source': 'phonopy',
                           'units':      {'frequency': units['frequency']},
                           'dimensions': {'frequency': dimensions['frequency']}}}

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

    for c in pconversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], float(pconversions[c]))

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], float(conversions[c]))

    return data2

def phonopy_gruneisen(filename):
    """Loads phonopy gruneisen data.

    Does not load path data, but can load from files with a q-point
    path, which will often be preferable if projecting onto a phonon
    dispersion.

    Arguments
    ---------

        filename : str
            phonopy gruneisen.yaml filepath.

    Returns
    -------

        dict
            gruneisen data.
    """

    import yaml

    pconversions = settings.phonopy_conversions()
    conversions = settings.conversions()

    # load data

    with open(filename, 'r') as f:
        data = yaml.safe_load(f)
    
    x = np.reshape([[q['distance'] for q in path['phonon']] for path in data['path']], (-1, 3))
    qp = np.reshape([[q['q-position'] for q in path['phonon']] for path in data['path']], (-1, 3))
    g = [[[band['gruneisen'] for band in q['band']] for q in path['phonon']] for path in data['path']]
    g = np.reshape(g, (-1, np.shape(g)[-1]))

    units = tp.settings.units()
    dimensions = settings.dimensions()
    data2 = {'x':             x,
             'qpoint':        qp,
             'gruneisen':     g,
             'meta':
                 {'gruneisen_source': 'phonopy',
                  'units':      {'gruneisen': units['gruneisen']},
                  'dimensions': {'gruneisen': dimensions['gruneisen']}}}

    for c in pconversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], float(pconversions[c]))

    for c in conversions:
        if c in data2:
            data2[c] = np.multiply(data2[c], float(conversions[c]))

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

def scale_to_path(x, tickpos, scalepos):
    """Scales data to a path.

    Useful to make different phonopy runs fit together or to map
    gruneisen data on a phonon dispersion.

    Arguments
    ---------

        x : list
            wavevector ordinates.
        tickpos : list
            tick wavevectors for scaling.
        scalepos : list
            scale tick wavevectors.

    Returns
    -------

        list
            wavevector ordinates.
    """

    n = 0
    # for each x, while within data and between two high-symmetry
    # points, interpolate onto scale data
    for i in range(len(x)):
        while n <= len(tickpos) and \
              not (x[i] >= tickpos[n] and x[i] <= tickpos[n+1]):
            n += 1
        x[i] = scalepos[n] + ((x[i] - tickpos[n]) * (scalepos[n+1] - scalepos[n]) / \
                             (tickpos[n+1] - tickpos[n]))

    return x
