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
import pymatgen as pmg
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Kpoints, Poscar

def gen_ibz(mesh, poscar='POSCAR', symprec=1e-5):
    """Generates the irreducible kpoints and weights for a mesh and material.

    Arguments
    ---------

        mesh : array-like
            3x1 kpoint mesh.

        poscar : str or pmg.io.vasp.inputs.Poscar, optional
            path to POSCAR file or pmg Poscar object. Default: POSCAR.
        symprec : float, optional
            symmetry precision. Default: 1e-5.

    Returns
    -------

        np.array
            irreducible kpoints.
        np.array
            weights.
    """

    if isinstance(poscar, str):
        sg = SpacegroupAnalyzer(Poscar.from_file(poscar).structure,
                                symprec=symprec)
    elif isinstance(poscar, pmg.io.vasp.inputs.Poscar):
        sg = SpacegroupAnalyzer(poscar.structure, symprec=symprec)
    ibz = sg.get_ir_reciprocal_mesh(mesh=mesh)

    kpts = np.array([list(k[0]) for k in ibz])
    weights = np.array([k[1] for k in ibz])

    return kpts, weights

def get_kpar(kpoints, poscar='POSCAR', symprec=1e-5):
    """Suggests KPARs for a set of kpoints and material.

    Ignores zero-weighted kpoints.

    Arguments
    ---------

        kpoints : array-like or str or pmg.io.vasp.inputs.Kpoints
            path to KPOINTS or IBZKPT file or 3x1 kpoint mesh or list of
            kpoint weights or list of kpoints or pmg Kpoints object.
            Lists of kpoints or weights should be 2D, and zero-weighted
            kpoints cannot be ignored in the case of lists of kpoints or
            KPOINTS files without weights.

        poscar : str or pmg.io.vasp.inputs.Poscar, optional
            path to POSCAR file or pmg Poscar object. Default: POSCAR.
        symprec : float, optional
            symmetry precision if mesh given. Default: 1e-5.

    Returns
    -------

        np.array
            potential KPARs in ascending order.
    """

    if isinstance(kpoints, (str, pmg.io.vasp.inputs.Kpoints)):
        if isinstance(kpoints, str): # KPOINTS file
            weights = Kpoints.from_file(kpoints).kpts_weights
            mesh = Kpoints.from_file(kpoints).kpts
        else: # pmg Kpoints object
            weights = kpoints.kpts_weights
            mesh = kpoints.kpts

        if weights is not None: # weighted kpoints
            weighted = len(np.array(weights).nonzero()[0])
        else:
            if np.shape(mesh) == (1, 3): # automatic mesh
                _, weights = gen_ibz(mesh, poscar, symprec)
                weighted = len(weights)

            else: # unweighted kpoints
                weighted = len(Kpoints.from_file(kpoints).kpts)

    elif len(np.shape(kpoints)) == 1: # mesh
        kpoints, weights = gen_ibz(kpoints, poscar, symprec)
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

    return [int(k) for k in np.sort(np.unique(kpar))]

def get_kpoints(kpoints, zero_weighted=None, poscar='POSCAR', output='KPOINTS'):
    """Generates a KPOINTS file.
    
       Includes optional zero-weighted kpoints. This will create a file
       if you supply it with KPOINTS meshes, however this is
       purportedly not necessarily the same as what vasp would do. To
       avoid this possibility, we recommend supplying two IBZKPT files,
       which this will stitch together.

    Arguments
    ---------

        kpoints : str or array-like or pmg.io.vasp.inputs.Kpoints
            IBZKPT file path or pmg Kpoints object or 3x1 kpoint mesh.

        zero-weighted : str or array-like or pmg.io.vasp.inputs.Kpoints, optional
            IBZKPT file path or pmg Kpoints object or 3x1 kpoint mesh.
            Default: no zero-weighted kpoints.
        poscar : str or pmg.io.vasp.inputs.Poscar, optional
            path to POSCAR file or pmg Poscar object. Default: POSCAR.
        output : str, optional
            output filename. Default: KPOINTS.

    Returns
    -------

        None
            writes to file
    """

    if isinstance(kpoints, (str, pmg.io.vasp.inputs.Kpoints)):
        known = False
        if isinstance(kpoints, str):
            f = Kpoints.from_file(kpoints)
        else:
            f = kpoints
        kpts = f.kpts
        weights = f.kpts_weights
    else:
        known = True
        kpts, weights = gen_ibz(kpoints, poscar)
    labels = list(np.full(len(kpts), ''))
    if known:
        labels[0] = '{} x {} x {} mesh'.format(kpoints[0], kpoints[1], kpoints[2])
    else:
        labels[0] = 'weighted mesh'

    if zero_weighted is not None:
        if isinstance(zero_weighted, (str, pmg.io.vasp.inputs.Kpoints)):
            known = False
            if isinstance(zero_weighted, str):
                zkpts = Kpoints.from_file(zero_weighted).kpts
            else:
                zkpts = zero_weighted.kpts
        else:
            known = True
            zkpts, _ = gen_ibz(zero_weighted, poscar)
        kpts = np.concatenate((kpts, zkpts))
        weights = np.concatenate((weights, np.zeros(len(zkpts))))
        zlabels = list(np.full(len(zkpts), ' '))
        if known:
            zlabels[0] = '{} x {} x {} zero-weighted mesh'.format(
                         zero_weighted[0], zero_weighted[1], zero_weighted[2])
        else:
            zlabels[0] = 'zero-weighted mesh'
        labels = np.concatenate((labels, zlabels))

    f = Kpoints(kpts=kpts, kpts_weights=list(weights), num_kpts=len(kpts),
                labels=list(labels), style='Reciprocal',
                comment='Generated with tp')
    f.write_file(output)

    return
