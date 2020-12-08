"""Utility to save to hdf5"""

import h5py
import numpy as np

def hdf5(data, output):
    """Saves to hdf5.

    Aims to make saving nested dictionaries easy, works for 3 layers.

    Arguments
    ---------

        data: dict
            data to save.
        output : str
            output filename.
    """

    datafile = h5py.File(output, 'w')

    for key in data.keys():
        if isinstance(data[key], dict):
            group = datafile.create_group(key)
            for k in data[key].keys():
                if isinstance(data[key][k], dict):
                    group2 = group.create_group(k)
                    for k2 in data[key][k].keys():
                        group2[k2] = data[key][k][k2]
                else:
                    group[k] = data[key][k]
        else:
            datafile.create_dataset(key, np.shape(data[key]), data=data[key])

    datafile.close()

    return
