"""Resolves quantities by temperature and/ or direction.

Reads variables and selects specific conditions. Requires
['meta']['dimensions'] subdictionaries provided by tp load modules.

Functions
---------

    merge:
        merge data from multiple files.
    resolve:
        selects data based on dependent variables.
"""

import numpy as np
import warnings
from copy import deepcopy

def merge(data, dependent):
    """Merges data dictionaries with shared dependent variables.

    Particularly with amset in mind, to relieve memory constraints.
    Currently works for one dependent variable, so looping will often
    bre required. An output in the form:

    ab

    cd

    would require three merges: a to b, c to d and ab to cd (or a to c etc.).

    Arguments
    ---------

        data : list of dicts
            data to merge. Requires tp metadata and dependent variable.
        dependent : string
            dependent variable.

    Returns
    -------

        dict
            merged data
    """
    depi = [list(range(len(data[0][dependent])))]
    data2 = {dependent: data[0][dependent],
             'meta': data[0]['meta']}
    for d in data[1:]:
        depi.append([d[dependent].index(d2) for d2 in d[dependent]\
                     if d2 not in data[0][dependent]])
        data2[dependent].extend(np.array(d[dependent])[depi[-1]])

    for key in data[0].keys():
        if key in [dependent, 'meta']:
            continue
        else:
            present = True
        for d in data[1:]:
            if key not in d.keys():
                present = False
                break
        if present and dependent in data[0]['meta']['dimensions'][key]:
            axis = data[0]['meta']['dimensions'][key].index(dependent)
            for d in data:
                d[key] = np.swapaxes(d[key], axis, 0)
            data2[key] = np.concatenate([np.array(data[i][key])[depi[i]] for i in range(len(data))], axis=0)
            for d in data:
                d[key] = np.swapaxes(d[key], 0, axis)
        elif key not in data2 and key in data[0]:
            data2[key] = data[0][key]
            
    return data2

def resolve(data, quantities, **kwargs):
    """Selects particular values of arbitrary quantities.

    Requires the meta/dimensions dictionaries found in later versions
    of tp. Currently cannot accept dictionary keys (e.g. dtype='n') if
    they are not in the 0th index.

    Arguments
    ---------

        data : dict
            data with meta/dimensions dictionaries and quantities.
        quantities : array-like or str
            quantities to resolve

        kwargs
            dimesions to resolve. Rounds to nearest available value.
            Common options include:

                direction
                    direction to resolve, accepts x-z/, a-c,
                    average/ avg/ mean/ arithmetic/ arith,  or
                    norm/ normal or harmonic/ harm.
                dtype
                    n or p
                stype
                    codes from amset, e.g. IMP, or overall
                doping
                    concentration, not to be confused with dtype
                temperature

    Returns
    -------

        dict
            resolved data.
    """

    data = deepcopy(data) # sever the link to enable the original data to be reused
    if 'meta' not in data or 'dimensions' not in data['meta']:
        raise Exception('data must have a meta subdictionary with a '
                        'dimensions subdictionary.')
    if isinstance(quantities, str):
        quantities = quantities.split()
    direction = {'a': 0, 'b': 1, 'c': 2,
                 'x': 0, 'y': 1, 'z': 2}

    # make sure dictionaries are dealt with first
    keys, vals = [], []
    for key, val in kwargs.items():
        if isinstance(val, str) and key != 'direction':
            keys.insert(0, key)
            vals.insert(0, val)
        else:
            keys.append(key)
            vals.append(val)

    for q in quantities:
        if q not in data['meta']['dimensions']:
            warnings.warn('{} not in dimensions. Skipping.'.format(q))
            continue
        for key, val in zip(keys, vals):
            if val is None:
                continue
            if key != 'direction':
                if key not in data and key not in ['dtype', 'stype']:
                    warnings.warn('{} not in data. Skipping.'.format(key))
                    continue
                if key not in data['meta']['dimensions'][q]:
                    continue
                if isinstance(val, str):
                    data['meta'][key] = val
                    for i, d in enumerate(data['meta']['dimensions'][q]):
                        if d == key:
                            if i == 0:
                                del data['meta']['dimensions'][q][i]
                                if key in data and val in data[key]:
                                    pos = data[key].index(val)
                                    data[q] = data[q][pos]
                                else:
                                    data[q] = data[q][val]
                                data['meta'][key] = val
                                break
                            else:
                                warnings.warn('Does not currently work '
                                              'unless strings are in '
                                              'the 0th index.')
                                break
                elif isinstance(val, (int, float, list, np.ndarray)):
                    if isinstance(val, (int, float)):
                        index = np.abs(np.subtract(data[key], val)).argmin()
                    else:
                        index = np.sqrt(np.sum(np.square(
                                np.subtract(data[key], val)), axis=1)).argmin()
                    data['meta'][key] = data[key][index]
                    for i, d in enumerate(data['meta']['dimensions'][q]):
                        if d == key:
                            del data['meta']['dimensions'][q][i]
                            data[q] = np.moveaxis(data[q], i, 0)
                            data[q] = data[q][index]
                            data['meta'][key] = data[key][index]
                            break
            else: # if key == 'direction':
                while True:
                    for i, d in enumerate(data['meta']['dimensions'][q]):
                        if d in [3, 6]:
                            del data['meta']['dimensions'][q][i]
                            data[q] = np.moveaxis(data[q], i, 0)
                            if val in direction:
                                data[q] = data[q][direction[val]]
                            elif val in ['mean', 'arithmetic', 'arith', 'average', 'avg']:
                                if len(data['meta']['dimensions'][q]) > i and \
                                   data['meta']['dimensions'][q][i] == 3:
                                    # if this is a 3x3 array
                                    del data['meta']['dimensions'][q][i]
                                    data[q] = np.moveaxis(data[q], i+1, 1)
                                    data[q] = np.average([data[q][0][0],
                                                          data[q][1][1],
                                                          data[q][2][2]],
                                                         axis=0)
                                else:
                                    # if this is a 3x1 or 6x1 array
                                    data[q] = np.average(data[q][:3], axis=0)
                                data['meta'][key] = 'arithmetic mean'
                            elif val in ['norm', 'normal']:
                                if len(data['meta']['dimensions'][q]) > i and \
                                   data['meta']['dimensions'][q][i] == 3:
                                    # if this is a 3x3 array
                                    del data['meta']['dimensions'][q][i]
                                    data[q] = np.moveaxis(data[q], i+1, 1)
                                    data[q] = np.square(data[q][0][0]) \
                                            + np.square(data[q][1][1]) \
                                            + np.square(data[q][2][2])
                                    data[q] = np.sqrt(data[q])
                                else:
                                    # if this is a 3x1 or 6x1 array
                                    data[q] = np.square(data[q][0]) \
                                            + np.square(data[q][1]) \
                                            + np.square(data[q][2])
                                    data[q] = np.sqrt(data[q])
                                data['meta'][key] = 'norm'
                            elif val in ['harmonic', 'harm']:
                                if len(data['meta']['dimensions'][q]) > i and \
                                   data['meta']['dimensions'][q][i] == 3:
                                    # if this is a 3x3 array
                                    del data['meta']['dimensions'][q][i]
                                    data[q] = np.moveaxis(data[q], i+1, 1)
                                    data[q] = 1/np.average([1/data[q][0][0],
                                                            1/data[q][1][1],
                                                            1/data[q][2][2]],
                                                           axis=0)
                                else:
                                    # if this is a 3x1 or 6x1 array
                                    data[q] = 1/np.average([1/data[q][0],
                                                            1/data[q][1],
                                                            1/data[q][2]],
                                                           axis=0)
                                data['meta'][key] = 'harmonic mean'
                            break
                    else:
                        break

    return data
