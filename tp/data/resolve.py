"""Resolves quantities by temperature and/ or direction.

Reads variables and selects specific conditions. Requires
['meta']['dimensions'] subdictionaries provided by tp load modules.

Functions
---------

    resolve
"""

import numpy as np
import tp.settings
import warnings
from copy import deepcopy

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
                    average/ avg or normal/ norm.
                dtype
                    n or p
                stype
                    codes from amset, e.g. IMP, or overall
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
                                data[q] = data[q][val]
                                data['meta'][key] = val
                                break
                            else:
                                warnings.warn('Does not currently work '
                                              'unless strings are in '
                                              'the 0th index.')
                                break
                elif isinstance(val, (int, float)):
                    index = np.abs(np.subtract(data[key], val)).argmin()
                    data['meta'][key] = data[key][index]
                    for i, d in enumerate(data['meta']['dimensions'][q]):
                        if d == key:
                            del data['meta']['dimensions'][q][i]
                            data[q] = np.moveaxis(data[q], i, 0)
                            data[q] = data[q][index]
                            data['meta'][key] = val
                            break
            else: # if key == 'direction'
                while True:
                    for i, d in enumerate(data['meta']['dimensions'][q]):
                        if d in [3, 6]:
                            del data['meta']['dimensions'][q][i]
                            data[q] = np.moveaxis(data[q], i, 0)
                            if val in direction:
                                data[q] = data[q][direction[val]]
                            elif val in ['avg', 'average']:
                                data[q] = np.average(data[q][:3], axis=0)
                            elif val in ['norm', 'normal']:
                                data[q] = np.square(data[q][0]) \
                                        + np.square(data[q][1]) \
                                        + np.square(data[q][2])
                                data[q] = np.sqrt(data[q])
                            data['meta'][key] = val
                            break
                    else:
                        break

    return data
