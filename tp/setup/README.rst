Contains functions to set up efficient VASP calculations.

Zero-weighted kpoints are an effective way to sample a larger kpoint
mesh at reduced computational cost: you need a converged weighted mesh,
and then a denser mesh  which is unweighted. ``get_kpoints`` will
combine two kpoint meshes, files or pymatgen objects to make a single
file with weighted and unweighted kpoints. Only weighted kpoints need
to be considered for KPAR, so ``get_kpar`` gives suggestions of KPAR
disregarding unweighted kpoints. ``get_ibz`` generates ibzkpoints to
help in these processes.
