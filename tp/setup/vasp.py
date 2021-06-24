"""Helps set up VASP calculations.

Functions
---------

    gen_ibz
        generate irreducible kpoints
    get_kpar
        gives kpar suggestions based on kpoints
    get_kpoints
        generates kpoints file including zero-weighted kpoints
"""

import numpy as np
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Kpoints, Poscar

def gen_ibz(mesh, poscar='POSCAR'):
    """Generates the irreducible kpoints and weights for a mesh and
       material.

    Arguments
    ---------

        mesh : array-like
            3x1 kpoint mesh.

        poscar : str, optional
            path to POSCAR file. Default: POSCAR.

    Returns
    -------

        np.array
            irreducible kpoints.
        np.array
            weights.
    """

    sg = SpacegroupAnalyzer(Poscar.from_file(poscar).structure)
    ibz = sg.get_ir_reciprocal_mesh(mesh=mesh)

    kpts = np.array([list(k[0]) for k in ibz])
    weights = np.array([k[1] for k in ibz])

    return kpts, weights

def get_kpar(kpoints, poscar='POSCAR'):
    """Suggests KPARs for a set of kpoints and material.

    Ignores zero-weighted kpoints.

    Arguments
    ---------

        kpoints : array-like or str
            path to KPOINTS or IBZKPT file or 3x1 kpoint mesh or list of
            kpoint weights or list of kpoints. Lists of kpoints or
            weights should be 2D, and zero-weighted kpoints cannot be
            ignored in the case of lists of kpoints or KPOINTS files
            without weights.

        poscar : str, optional
            path to POSCAR file. Default: POSCAR.

    Returns
    -------

        np.array
            potential KPARs in ascending order.
    """

    if isinstance(kpoints, str): # read from file
        weights = Kpoints.from_file(kpoints).kpts_weights
        if weights is not None: # weighted kpoints
            weighted = len(np.array(weights).nonzero()[0])
        else:
            mesh = Kpoints.from_file(kpoints).kpts
            if np.shape(mesh) == (1, 3): # automatic mesh
                sg = SpacegroupAnalyzer(Poscar.from_file(poscar).structure)
                ibz = sg.get_ir_reciprocal_mesh(mesh=mesh)
                weighted = len(ibz)
            else: # unweighted kpoints
                weighted = len(Kpoints.from_file(kpoints).kpts)
    elif len(np.shape(kpoints)) == 1: # mesh
        kpoints, weights = gen_ibz(kpoints, poscar)
        weighted = len(weights.nonzero()[0])
    elif len(np.shape(kpoints)) == 2: # list...
        if len(kpoints[0]) == 3: # ...of kpoints
            weighted = len(kpoints)
        else: # ...of weights
            weighted = len(np.array(kpoints).nonzero()[0])

    # finds the factors of the number of kpoints
    kpar = []
    for n in range(1, int(np.ceil(np.sqrt(weighted))) + 1):
        if weighted % n == 0:
            kpar.append(n)
            kpar.append(weighted / n)

    return np.sort(np.unique(kpar))

def get_kpoints(mesh, zero_weighted=None, poscar='POSCAR', output='KPOINTS'):
    """Generates a KPOINTS file, including optional zero-weighted
       kpoints, for kpoint mesh(es) and a material.

    Arguments
    ---------

        kpoints : array-like
            3x1 kpoint mesh.

        zero-weighted : array-like, optional
            3x1 kpoint mesh. Default: no zero-weighted kpoints.
        poscar : str, optional
            path to POSCAR file. Default: POSCAR.
        output : str, optional
            output filename. Default: KPOINTS.

    Returns
    -------

        None
            writes to file
    """

    kpts, weights = gen_ibz(mesh, poscar)
    labels = list(np.full(len(kpts), ' '))
    labels[0] = '{} x {} x {} mesh'.format(mesh[0], mesh[1], mesh[2])

    if zero_weighted is not None:
        zkpts, _ = gen_ibz(zero_weighted, poscar)
        kpts = np.concatenate((kpts, zkpts))
        weights = np.concatenate((weights, np.zeros(len(zkpts))))
        zlabels = list(np.full(len(zkpts), ' '))
        zlabels[0] = '{} x {} x {} zero-weighted mesh'.format(
                     zero_weighted[0], zero_weighted[1], zero_weighted[2])
        labels = np.concatenate((labels, zlabels))

    f = Kpoints(kpts=kpts, kpts_weights=list(weights), num_kpts=len(kpts),
                labels=list(labels), style='Reciprocal',
                comment='Generated with tp')
    f.write_file(output)

    return
