"""Resolves quantities by temperature and/ or q-point index.

May need to be split by data origin in future.

Functions:
    resolve:
        currently for Phono3py, parts of AMSET and BoltzTraP.
"""

import numpy as np
import tp.settings
from tp.data import aniso

def resolve(data, quantities, temperature=None, direction=None):
    """Selects temperature and/or direction.

    AMSET properties currently only support direction.

    Arguments:
        data : dict
            data.
        quantities : array-like or str
            quantities to resolve

        temperature : float, optional
            temperature to select. Default: None.
        direction : str, optional
            direction to resolve, accepts x-z/, a-c, average/ avg or
            normal/ norm. Default: None.

    Returns:
        dict
            resolved data.
    """

    data = dict(data) # sever the link to enable the original data to be reused

    # variables resolved by temperature and direction
    hast = ['average_eff_mass', 'conductivity',
            'electronic_thermal_conductivity', 'fermi_level', 'gamma',
            'heat_capacity', 'lattice_thermal_conductivity', 'lifetime',
            'mean_free_path', 'mode_kappa', 'occupation', 'power_factor',
            'seebeck', 'velocities_product', 'zt']
    iso = {'average_eff_mass':                aniso.matrix_three,
           'conductivity':                    aniso.matrix_three,
           'electronic_thermal_conductivity': aniso.matrix_three,
           'group_velocity':                  aniso.three,
           'gv_by_gv':                        aniso.three,
           'lattice_thermal_conductivity':    aniso.two,
           'mean_free_path':                  aniso.four,
           'mesh':                            aniso.one,
           'mode_kappa':                      aniso.four,
           'power_factor':                    aniso.matrix_three,
           'qpoint':                          aniso.two,
           'seebeck':                         aniso.matrix_three,
           'velocities_product':              aniso.matrix_two,
           'zt':                              aniso.matrix_three}

    # temperature resolution

    if temperature is not None and 'temperature' in data:
        ti = np.abs(np.array(data['temperature']) - temperature).argmin()
        data['meta']['temperature'] = data['temperature'][ti]

        iso['average_eff_mass'] = aniso.matrix_two
        iso['conductivity'] = aniso.matrix_two
        iso['electronic_thermal_conductivity'] = aniso.matrix_two
        iso['lattice_thermal_conductivity'] = aniso.one
        iso['mean_free_path'] = aniso.three
        iso['mode_kappa'] = aniso.three
        iso['power_factor'] = aniso.matrix_two
        iso['seebeck'] = aniso.matrix_two
        iso['velocities_product'] = aniso.matrix_one
        iso['zt'] = aniso.matrix_two
    if direction is not None:
        data['meta']['direction'] = direction

    # direction resolution

    tnames = tp.settings.to_tp()
    if isinstance(quantities, str):
        quantities = quantities.split()
    quantities2 = []
    for i in quantities:
        if i not in quantities2:
            quantities2.append(i)
    for q in quantities2:
        q2 = tnames[q] if q in tnames else q
        if temperature is not None and q2 in hast:
            data[q] = data[q][ti]
        if direction is not None and q2 in iso:
            data[q] = iso[q2](data[q], direction)

    return data
