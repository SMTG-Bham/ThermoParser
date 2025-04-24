"""Helps set up phonopy calculations."""

#Functions
#---------
#
#    gen_band_conf
#        generate phonopy band.conf.
#    gen_dos_conf
#        generate phonopy dos.conf.
#"""

import numpy as np
import pymatgen as pmg
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.kpath import KPathSeek, KPathLatimerMunro, KPathSetyawanCurtarolo
import re

def gen_band_conf(dim, preset='pymatgen', path=None, labels=None,
                  symprec=1e-5, points=101, connection=True,
                  eigenvectors=True, nac=False, poscar='POSCAR',
                  output='band.conf'):
    """Generates a band.conf file to run phonopy with.

    Arguments
    ---------

        dim : array-like
            isotropic scale factor or 3x1 or 3x3 or 6x1 supercell size matrix.
        preset : str, optional
            path formalism. Choice of seekpath/hinuma, pymatgen/setyawan,
            latimer and bradley. bradley requires sumo to be installed.
            Overridden by path. Default: pymatgen
        path : str or array-like, optional
            custom k-path. Overrides preset.
        labels : str or arraylike, optional
            custom label names for custom q-paths.
        symprec : float, optional
            symmetry tolerance for space group identification in AA.
            Default: 1e-5
        points : int, optional
            number of points per path. Default: 101
        connection : bool, optional
            attempt to maintain band order at high-symmetry points.
            Default: True
        eigenvectors : bool, optional
            print eigenvalues to file. Default: True
        nac : bool, optional
            apply non-analytical correction. Requires phonopy BORN file.
            Default: False

        poscar : str or pmg.io.vasp.inputs.Poscar, optional
            path to POSCAR file or pmg Poscar object. Default: POSCAR
        output : str, optional
            output name. Default: band.conf

    Returns
    -------

        None
            prints to output
    """

    if isinstance(dim, str):
        dim = dim.split()
    if np.shape(dim) == (3,) or np.shape(dim) == (9,):
        dim = [str(d) for d in dim]
        dimstr = ' '.join(dim)
    elif np.shape(dim) == (6,):
        dimstr = f'{dim[0]} {dim[5]} {dim[4]} {dim[5]} {dim[1]} {dim[3]} {dim[4]} {dim[3]} {dim[2]}'
    elif np.shape(dim) == (3,3):
        dimstr = f'{dim[0][0]} {dim[0][1]} {dim[0][2]} {dim[1][0]} {dim[1][1]} {dim[1][2]} {dim[2][0]} {dim[2][1]} {dim[2][2]}'
    elif np.shape(dim) == (1,):
        dimstr = f'{dim[0]} {dim[0]} {dim[0]}'
    else:
        raise ValueError("dim must be a 1x1 or 3x1 or 3x3 or 6x1 or 9x1 array")

    if path is None:
        ps = {'bradley_cracknell': 'bradley',
              'bradley':           'bradley',
              'bradcrack':         'bradley',
              'hinuma':            'hinuma',
              'seekpath':          'hinuma',
              'latimer_munro':     'latimer',
              'latimer':           'latimer',
              'setyawan_cutarolo': 'setyawan',
              'setyawan':          'setyawan',
              'pymatgen':          'setyawan'
              }
        if isinstance(poscar, str):
            poscar = Poscar.from_file(poscar).structure
        elif isinstance(poscar, pmg.io.vasp.inputs.Poscar):
            poscar = poscar.structure

        if path is None:
            if ps[preset] == 'bradley':
                try:
                    from sumo.symmetry.brad_crack_kpath import BradCrackKpath
                except ModuleNotFoundError:
                    raise Exception("sumo is required for Bradley-Cracknell kpaths. "
                                    "e.g.: pip install sumo")
                kpath = BradCrackKpath(poscar, symprec=symprec)._kpath
            else:
                psfunc = {'hinuma':   KPathSeek,
                          'latimer':  KPathLatimerMunro,
                          'setyawan': KPathSetyawanCurtarolo}
                kpath = psfunc[ps[preset]](poscar, symprec=symprec)._kpath
        else:
            kpath = {}
            if isinstance(path, str):
                path = ' '.split(path)
            kpath['kpoints'] = np.reshape(kpath, (-1, 3))
            if labels is not None:
                if isinstance(labels, str):
                    kpath['path'] = ' '.split(labels)
                else:
                    kpath['path'] = labels
        
        if 'path' in kpath:
            bandstr = ''
            labelstr = ''
            for i, path in enumerate(kpath['path']):
                for j, point in enumerate(path):
                    labelpoint = point + ' '
                    bandpoint = str(kpath['kpoints'][point][0]) + ' '
                    bandpoint += str(kpath['kpoints'][point][1]) + ' '
                    bandpoint += str(kpath['kpoints'][point][2]) + '  '
                    labelpoint = re.sub(r'\\Gamma', r'$\\Gamma$', labelpoint)
                    labelpoint = re.sub('GAMMA', r'$\\Gamma$', labelpoint)
                    labelpoint = re.sub('Γ', r'$\\Gamma$', labelpoint)
                    while len(labelpoint) < len(bandpoint):
                        labelpoint += ' ' 
                    while len(labelpoint) > len(bandpoint):
                        bandpoint += ' ' 
                    bandstr += bandpoint
                    labelstr += labelpoint
                if i+1 < len(kpath['path']):
                    bandstr = re.sub(' * $', ' , ', bandstr)
                    #nextpoint = kpath['path'][i+1][0]
                    #labelstr = re.sub(f'{point} *$', f'"{point} | {nextpoint}" ', labelstr)
                    while len(labelstr) < len(bandstr):
                        labelstr += ' ' 
                    while len(labelstr) > len(bandstr):
                        bandstr += ' ' 

    with open(output, 'w') as f:
        f.write(f'DIM = {dimstr}\n')
        f.write(f'BAND        = {bandstr}\n')
        if 'path' in kpath:
            f.write(f'BAND_LABELS = {labelstr}\n')
        f.write(f'BAND_POINTS = {points}\n')
        f.write(f'BAND_CONNECTION = .{str(connection).upper()}.\n')
        f.write(f'EIGENVECTORS = .{str(eigenvectors).upper()}.\n')
        f.write(f'NAC = .{str(nac).upper()}.')

    return kpath

def gen_dos_conf(dim, mesh=[24, 24, 24], eigenvectors=True,
                 gamma=True, nac=False, poscar='POSCAR',
                 output='dos.conf'):
    """Generates a dos.conf file to run phonopy with.

    Arguments
    ---------

        dim : array-like
            3x1 or 3x3 supercell size.
        mesh : array-like, optional
            3x1 qpoint mesh. Default: [24, 24, 24]
        eigenvectors : bool, optional
            print eigenvalues to file. Default: True
        gamma : bool, optional
            use gamma-centred q-mesh. Default: True
        nac : bool, optional
            apply non-analytical correction. Requires phonopy BORN file.
            Default: False

        poscar : str or pmg.io.vasp.inputs.Poscar, optional
            path to POSCAR file or pmg Poscar object. Default: POSCAR
        output : str, optional
            output name. Default: dos.conf

    Returns
    -------

        None
            prints to output
    """

    if isinstance(dim, str):
        dim = dim.split()
    if np.shape(dim) == (3,) or np.shape(dim) == (9,):
        dim = [str(d) for d in dim]
        dimstr = ' '.join(dim)
    elif np.shape(dim) == (6,):
        dimstr = f'{dim[0]} {dim[5]} {dim[4]} {dim[5]} {dim[1]} {dim[3]} {dim[4]} {dim[3]} {dim[2]}'
    elif np.shape(dim) == (3,3):
        dimstr = f'{dim[0][0]} {dim[0][1]} {dim[0][2]} {dim[1][0]} {dim[1][1]} {dim[1][2]} {dim[2][0]} {dim[2][1]} {dim[2][2]}'
    elif np.shape(dim) == (1,):
        dimstr = f'{dim[0]} {dim[0]} {dim[0]}'
    else:
        raise ValueError("dim must be a 1x1 or 3x1 or 3x3 or 6x1 or 9x1 array")

    if isinstance(mesh, str):
        mesh = mesh.split()
    if np.shape(mesh) == (3,):
        mesh = [str(m) for m in mesh]
        meshstr = ' '.join(mesh)
    else:
        raise ValueError("mesh must be a 3x1 array")

    if isinstance(poscar, str):
        poscar = Poscar.from_file(poscar).structure
    elif isinstance(poscar, pmg.io.vasp.inputs.Poscar):
        poscar = poscar.structure

    atomlist = ''
    xatom = None
    numberlist = ''
    for i, atom in enumerate(poscar):
        if atom.specie != xatom:
            if xatom is not None:
                numberlist += ', '
                while len(numberlist) < len(atomlist):
                    numberlist += ' ' 
                while len(numberlist) > len(atomlist):
                    atomlist += ' ' 
            xatom = atom.specie
            atomlist += str(atom.specie) + ' '
        numberlist += str(i+1) + ' '

    with open(output, 'w') as f:
        f.write(f'DIM = {dimstr}\n')
        f.write(f'MP = {meshstr}\n')
        f.write(f'PDOS  = {numberlist}\n')
        f.write(f'ATOMS = {atomlist}\n')
        f.write(f'GAMMA_CENTER = .{str(gamma).upper()}.\n')
        f.write(f'EIGENVECTORS = .{str(eigenvectors).upper()}.\n')
        f.write(f'NAC = .{str(nac).upper()}.')

    return