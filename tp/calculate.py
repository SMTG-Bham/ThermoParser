"""Calculates stuff.

Uses tp units by default, but can read the tprc.yaml to customise.

Functions
---------

    cumulate:
        sorts and cumulates.
    lorentzian:
        lorentzian curve.
    lifetime:
        particle lifetime.
    mfp:
        particle mean free path.
    be_occupation:
        boson occupation.
    dfdde:
        derivative of the Fermi-Dirac distribution.
    power_factor:
        power factor.
    zt:
        ZT.
    kl:
        lattice thermal conductivity for target ZT.

    power_factor_fromdict:
        adds power factor to dictionary.
    zt_fromdict:
        adds zt to dictionary.
    kl_fromdict:
        adds lattice thermal conductivity for target ZT to dictionary.

    to_tp
        converts quantities to tp defaults from tprc.yaml.
    from_tp
        converts quantities from tp defaults from tprc.yaml.
"""

import numpy as np
import tp
from scipy.constants import physical_constants

def cumulate(x, y):
    """Sorts by x and cumulates y.

    Arguments
    ---------

        x : array-like
            x-values.
        y : array-like
            y-values.

    Returns
    -------

        np.array
            sorted x-values.
        np.array
            cumulated y-values.
    """

    x = np.ravel(x)
    xsort = x[x.argsort()]

    y = np.ravel(y)
    ysort = y[x.argsort()]
    ycum = np.cumsum(np.ma.masked_invalid(ysort))

    return xsort, ycum

def lorentzian(x, x0=0, fwhm=1):
    """Area conserved Lorentzian function centered on x0.

    Arguments
    ---------

        x : array-like
            x-values.
        x0 : float
            origin of function.
        fwhm : float
            full-width at half-maximum.

    Returns
    -------

        np.array
            lorentzian
    """

    x = np.array(x)

    return 0.5 * fwhm / (np.pi * ((x - x0)**2 + (0.5 * fwhm)**2))

def lifetime(gamma, use_tprc=True):
    """Calculates lifetime from imaginary self-energy (Gamma).

    Arguments
    ---------

        gamma : array-like or float
            frequencies (by default in THz).

        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        np.array
            lifetimes (by default in s).
    """

    if use_tprc:
        gamma = to_tp('gamma', gamma)

    lifetime = np.reciprocal(np.multiply(2 * 2 * np.pi, gamma))
    lifetime = np.where(np.isinf(lifetime), 0, data['lifetime'])

    if use_tprc:
        lifetime = from_tp('lifetime', lifetime)

    return lifetime

def mfp(gamma, group_velocity, use_tprc=True):
    """
    Calculates mean free path from imaginary self-energy (Gamma) and
    group velocity.

    Arguments
    ---------

        gamma : array-like or float
            frequencies (by default in THz).
        group_velocity : array-like or float
            group velocities (by default in m s-1).

        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        array-like
            mean free paths (by default in m).
    """

    if use_tprc:
        group_velocity = to_tp('group_velocity', group_velocity)

    tau = lifetime(gamma, use_tprc=use_tprc)
    mfp = np.multiply(np.transpose([tau,] * 3, (1,2,3,0)), group_velocity)

    if use_tprc:
        mfp = from_tp('mean_free_path', mfp)

    return mfp

def be_occupation(frequency, temperature, use_tprc=True):
    """Calculates Bose-Einstein occupation.

    Arguments
    ---------

        frequency : array-like or float
            frequencies (by default in THz).
        temperature : array-like or float
            temperature (by default in K).

        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        array-like
            occupations.
    """

    hbar = physical_constants['Planck constant over 2 pi in eV s'][0]
    kb = physical_constants['Boltzmann constant in eV/K'][0]

    if use_tprc:
        frequency = to_tp('frequency', frequency)
        temperature = to_tp('temperature', temperature)

    frequency = np.array(frequency)
    occupation = np.expm1(np.divide.outer(kb * temperature,
                                          frequency * 1e12 * hbar) ** -1) ** -1

    if use_tprc:
        occupation = from_tp('occupation', occupation)

    return occupation

def dfdde(energy, fermi_levels, temperature, doping, amset_order=False,
          use_tprc=True):
    """Derivative of the Fermi-Dirac distribution.

    Arguments
    ---------

        energy : array-like
            energies per band and k-point (by default in eV).
        fermi_levels : array-like
            fermi levels per temperature and dopant (by default in eV).
        temperature : array-like
            temperatures (by default in K).
        doping : array-like
            doping concentrations (by default in cm-3).

        amset_order : bool, optional
            doping index before temperature index. Default: False.
        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        np.array
            derivative of the Fermi-Dirac distribution.
    """

    kb = physical_constants['Boltzmann constant in eV/K'][0]

    if use_tprc:
        energy = to_tp('energy', energy)
        fermi_level = to_tp('energy', fermi_level)
        temperature = to_tp('temperature', temperature)
        doping = to_tp('doping', doping)

    kbt = np.multiply(kb, temperature)
    de = -np.subtract.outer(fermi_levels, energy)
    if amset_order:
        weights = -0.25 / np.cosh(0.5 * de / kbt[None, :, None, None]) ** 2
        weights = weights / kbt[None, :, None, None]
    else:
        weights = -0.25 / np.cosh(0.5 * de / kbt[:, None, None, None]) ** 2
        weights = weights / kbt[:, None, None, None]

    return weights

def power_factor(conductivity, seebeck, use_tprc=True):
    """Calculates power factor.

    Arguments
    ---------

        conductivity : array-like
            conductivities (by default in S m-1).
        seebeck : array-like
            seebeck coefficients (by default in muV K-1).

        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        np.array
            power factors (by default in W m-1 K-2).
    """

    if use_tprc:
        conductivity = to_tp('conductivity', conductivity)
        seebeck = to_tp('seebeck', seebeck)

    pf = np.multiply(conductivity, 1e-12 * np.square(seebeck))

    if use_tprc:
        pf = from_tp('power_factor', pf)

    return pf

def zt(conductivity, seebeck, electronic_thermal_conductivity,
       lattice_thermal_conductivity, temperature, use_tprc=True):
    """Calculates ZT.

    Arguments
    ---------

        conductivity : array-like
            conductivities (by default in S m-1).
        seebeck : array-like
            seebeck coefficients (by default in muV K-1).
        electronic_thermal_conductivity : array-like
            electronic thermal conductivities (by default in W m-1 K-1).
        lattice_thermal_conductivity : array-like
            lattice thermal conductivities (by default in W m-1 K-1).
        temperature : array-like
            temperatures (by default in K).

        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        np.array
            ZT.
    """

    pf = power_factor(conductivity, seebeck, use_tprc=use_tprc)

    if use_tprc:
        pf = to_tp('power_factor', pf)
        electronic_thermal_conductivity = to_tp(
            'electronic_thermal_conductivity', electronic_thermal_conductivity)
        lattice_thermal_conductivity = to_tp('lattice_thermal_conductivity',
                                              lattice_thermal_conductivity)
        temperature = to_tp('temperature', temperature)

    zt = np.multiply(pf, np.array(temperature)[:, None]) / \
                     np.add(electronic_thermal_conductivity,
                            np.array(lattice_thermal_conductivity)[:, None])

    if use_tprc:
        zt = from_tp('zt', zt)

    return zt

def kl(conductivity, seebeck, electronic_thermal_conductivity, zt, temperature,
       use_tprc=True):
    """Calculates lattice thermal conductivity.

    Arguments
    ---------

        conductivity : array-like
            conductivities (by default in S m-1).
        seebeck : array-like
            seebeck coefficients (by default in muV K-1).
        electronic_thermal_conductivity : array-like
            electronic thermal conductivities (by default in W m-1 K-1).
        zt : array-like
            zt.
        temperature : array-like
            temperatures (by default in K).

        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        np.array
            lattice thermal conductivities (by default in W m-1 K-1).
    """

    pf = power_factor(conductivity, seebeck, use_tprc=use_tprc)

    if use_tprc:
        pf = to_tp('power_factor', pf)
        electronic_thermal_conductivity = to_tp(
            'electronic_thermal_conductivity', electronic_thermal_conductivity)
        zt = to_tp('zt', zt)
        temperature = to_tp('temperature', temperature)

    mid = np.divide(pf * np.array(temperature)[:, None], zt)
    kl = np.subtract(mid, electronic_thermal_conductivity)

    if use_tprc:
        kl = from_tp('lattice_thermal_conductivity', kl)

    return kl

def power_factor_fromdict(data, use_tprc=True):
    """Convenience wrapper to calculate power factor from a dictionary.

    Arguments
    ---------

        data : dict
            dictionary containing:

                conductivity : array-like
                    conductivities (by default in S m-1).
                seebeck : array-like
                    seebeck coefficients (by default in muV K-1).

        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        dict
            dictionary with power factors (by default in W m-1 K-2).
    """

    data['power_factor'] = power_factor(data['conductivity'],
                                        data['seebeck'], use_tprc=use_tprc)
    data['meta']['units']['power_factor'] = \
                           tp.settings.units(use_tprc=use_tprc)['power_factor']

    return data

def zt_fromdict(data, use_tprc=True):
    """Convenience wrapper to calculate ZT from a dictionary.

    Arguments
    ---------

        data : dict
            dictionary containing:

                conductivity : array-like
                    conductivities (by default in S m-1).
                seebeck : array-like
                    seebeck coefficients (by default in muV K-1).
                electronic_thermal_conductivity : array-like
                    electronic thermal conductivities (by default
                    in W m-1 K-1).
                lattice_thermal_conductivity : array-like
                    lattice thermal conductivities (by default
                    in W m-1 K-1).
                temperature : array-like
                    temperatures (by default in K).

        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        dict
            dictionary with ZTs.
    """

    data['zt'] = zt(data['conductivity'], data['seebeck'],
                 data['electronic_thermal_conductivity'],
                 data['lattice_thermal_conductivity'],
                 data['temperature'], use_tprc=use_tprc)
    data['meta']['units']['zt'] = tp.settings.units(use_tprc=use_tprc)['zt']

    return data

def kl_fromdict(data, use_tprc=True):
    """Convenience wrapper to calculate k_latt from a dictionary.

    Arguments
    ---------

        data : dict
            dictionary containing:

                conductivity : array-like
                    conductivities (by default in S m-1).
                seebeck : array-like
                    seebeck coefficients (by default in muV K-1).
                electronic_thermal_conductivity : array-like
                    electronic thermal conductivities (by default
                    in W m-1 K-1).
                zt : array-like
                    zt.
                temperature : array-like
                    temperatures (by default in K).

        use_tprc : bool, optional
            use custom unit conversions. Default: True.

    Returns
    -------

        dict
            dictionary with lattice thermal conductivities (by default
            in W m-1 K-1).
    """

    q = 'lattice_thermal_conductivity'
    data[q] = kl(data['conductivity'], data['seebeck'],
                 data['electronic_thermal_conductivity'], data['zt'],
                 data['temperature'], use_tprc=use_tprc)
    data['meta']['units'][q] = tp.settings.units(use_tprc=use_tprc)[q]

    return data

def to_tp(name, value):
    """Converts quantity to tp default units using tprc.yaml.

    Arguments
    ---------

        value : array-like
            values to be converted.
        name : str
            name in tprc.yaml.

    Returns
    -------

        np.array
            converted values.
    """

    conversions = tp.settings.conversions()
    if name in conversions and conversions[name] is not None:
        value = np.divide(value, conversions[name])

    return value

def from_tp(name, value):
    """Converts quantity from tp default units using tprc.yaml.

    Arguments
    ---------

        value : array-like
            values to be converted.
        name : str
            name in tprc.yaml.

    Returns
    -------

        np.array
            converted values.
    """

    conversions = tp.settings.conversions()
    if name in conversions and conversions[name] is not None:
        value = np.multiply(value, conversions[name])

    return value
