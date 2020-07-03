"""Functions for resolving anisotropic variables.

Naming convention: number indicates ordinal to be resolved; a lone
number indicates the form [x, y, z] (or [xx, yy, zz, yz, xz, xy] etc.),
while matrix_number indicates a three by three matrix.
"""

import numpy as np

def one(data, direction='avg'):
    """Resolves first-index anisotropy.

    Arguments:
        data : array-like
            anisotropic variable.
        direction : str, optional
            direction to resolve. Accepts a-c/ x-y, average/ avg. or
            normal/ norm. Default: average.

    Returns:
        list
            resolved variable.
    """

    def a(data): return np.array(data)[0]
    def b(data): return np.array(data)[1]
    def c(data): return np.array(data)[2]
    def avg(data): return np.mean([a(data), b(data), c(data)], axis=0)
    def norm(data):
        n = np.sqrt(np.sum(np.power([a(data), b(data), c(data)], 2), axis=0))
        return n

    dirn = {'a': a, 'b': b, 'c': c, 'x': a, 'y': b, 'z': c,
            'avg': avg, 'average': avg, 'norm': norm, 'normal': norm}

    assert direction in dirn, '{} unrecognised. Direction must be in {}.'.format(
                              direction, ', '.join(dirn))

    data2 = dirn[direction](data)

    return data2

def two(data, direction='avg'):
    """Resolves second-index anisotropy.

    Arguments:
        data : array-like
            anisotropic variable.
        direction : str, optional
            direction to resolve. Accepts a-c/ x-y, average/ avg. or
            normal/ norm. Default: average.

    Returns:
        list
            resolved variable.
    """

    def a(data): return np.array(data)[:,0]
    def b(data): return np.array(data)[:,1]
    def c(data): return np.array(data)[:,2]
    def avg(data): return np.mean([a(data), b(data), c(data)], axis=0)
    def norm(data):
        n = np.sqrt(np.sum(np.power([a(data), b(data), c(data)], 2), axis=0))
        return n

    dirn = {'a': a, 'b': b, 'c': c, 'x': a, 'y': b, 'z': c,
            'avg': avg, 'average': avg, 'norm': norm, 'normal': norm}

    assert direction in dirn, '{} unrecognised. Direction must be in {}.'.format(
                              direction, ', '.join(dirn))

    return dirn[direction](data)

def three(data, direction='avg'):
    """Resolves third-index anisotropy.

    Arguments:
        data : array-like
            anisotropic variable.
        direction : str, optional
            direction to resolve. Accepts a-c/ x-y, average/ avg. or
            normal/ norm. Default: average.

    Returns:
        list
            resolved variable.
    """

    def a(data): return np.array(data)[:,:,0]
    def b(data): return np.array(data)[:,:,1]
    def c(data): return np.array(data)[:,:,2]
    def avg(data): return np.mean([a(data), b(data), c(data)], axis=0)
    def norm(data):
        n = np.sqrt(np.sum(np.power([a(data), b(data), c(data)], 2), axis=0))
        return n

    dirn = {'a': a, 'b': b, 'c': c, 'x': a, 'y': b, 'z': c,
            'avg': avg, 'average': avg, 'norm': norm, 'normal': norm}

    assert direction in dirn, '{} unrecognised. Direction must be in {}.'.format(
                              direction, ', '.join(dirn))

    return dirn[direction](data)

def four(data, direction='avg'):
    """Resolves fourth-index anisotropy.

    Arguments:
        data : array-like
            anisotropic variable.
        direction : str, optional
            direction to resolve. Accepts a-c/ x-y, average/ avg. or
            normal/ norm. Default: average.

    Returns:
        list
            resolved variable.
    """

    def a(data): return np.array(data)[:,:,:,0]
    def b(data): return np.array(data)[:,:,:,1]
    def c(data): return np.array(data)[:,:,:,2]
    def avg(data): return np.mean([a(data), b(data), c(data)], axis=0)
    def norm(data):
        n = np.sqrt(np.sum(np.power([a(data), b(data), c(data)], 2), axis=0))
        return n

    dirn = {'a': a, 'b': b, 'c': c, 'x': a, 'y': b, 'z': c,
            'avg': avg, 'average': avg, 'norm': norm, 'normal': norm}

    assert direction in dirn, '{} unrecognised. Direction must be in {}.'.format(
                              direction, ', '.join(dirn))

    return dirn[direction](data)

def matrix_three(data, direction='avg'):
    """Resolves third-index anisotropy for matrix representations.

    Arguments:
        data : array-like
            anisotropic variable.
        direction : str, optional
            direction to resolve. Accepts a-c/ x-y, average/ avg. or
            normal/ norm. Default: average.

    Returns:
        list
            resolved variable.
    """

    def a(data): return np.array(data)[:,:,0,0]
    def b(data): return np.array(data)[:,:,1,1]
    def c(data): return np.array(data)[:,:,2,2]
    def avg(data): return np.mean([a(data), b(data), c(data)], axis=0)
    def norm(data):
        n = np.sqrt(np.sum(np.power([a(data), b(data), c(data)], 2), axis=0))
        return n

    dirn = {'a': a, 'b': b, 'c': c, 'x': a, 'y': b, 'z': c,
            'avg': avg, 'average': avg, 'norm': norm, 'normal': norm}

    assert direction in dirn, '{} unrecognised. Direction must be in {}.'.format(
                              direction, ', '.join(dirn))

    return dirn[direction](data)
