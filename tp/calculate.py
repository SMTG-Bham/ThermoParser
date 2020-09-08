"""Calculates stuff.

Functions:
    cumulate:
        sorts and cumulates.
    lorentzian:
        lorentzian curve.
    be_occupation:
        boson occupation.
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
"""

import numpy as np
import tp

def cumulate(x, y):
    """Sorts by x and cumulates y."""

    x = np.ravel(x)
    xsort = x[x.argsort()]

    y = np.ravel(y)
    ysort = y[x.argsort()]
    ycum = np.cumsum(np.ma.masked_invalid(ysort))

    return xsort, ycum

def lorentzian(x, x0=0, fwhm=1):
    """Area conserved Lorentzian function centered on x0.

    Arguments:
        x : np.array
            x-values.
        x0 : float
            origin of function.
        fwhm : float
            full-width at half-maximum.

    Returns:
        np.array
            lorentzian
    """

    x = np.array(x)

    return 0.5 * fwhm / (np.pi * ((x - x0)**2 + (0.5 * fwhm)**2))

def be_occupation(frequency, temperature=300.):
    """Calculates Bose-Einstein occupation.

    Arguments:
        frequency : array-like or float
            frequencies in THz.

        temperature : float, optional
            temperature in K. Default: 300.

    Returns:
        array-like
            occupations.
    """

    import scipy.constants as con

    hbar = con.physical_constants['Planck constant over 2 pi in eV s'][0]
    kb = con.physical_constants['Boltzmann constant in eV/K'][0]

    frequency = np.array(frequency)
    occupation = np.expm1((frequency * 1e12 * hbar) / (kb * temperature)) ** -1

    return occupation

def power_factor(conductivity, seebeck):
    """Calculates power factor.

    Arguments:
        conductivity : array-like
            conductivities.
        seebeck : array-like
            seebeck coefficients.

    Returns:
        np.array
            power factors.
    """

    return np.multiply(conductivity, 1e-12 * np.square(seebeck))

def zt(conductivity, seebeck, electronic_thermal_conductivity,
       lattice_thermal_conductivity, temperature):
    """Calculates ZT.

    Arguments:
        conductivity : array-like
            conductivities.
        seebeck : array-like
            seebeck coefficients.
        electronic_thermal_conductivity : array-like
            electronic thermal conductivities.
        lattice_thermal_conductivity : array-like
            lattice thermal conductivities by temperature.
        temperature : array-like
            temperatures.

    Returns:
        np.array
            ZT.
    """

    pf = power_factor(conductivity, seebeck)
    zt = np.multiply(pf, np.array(temperature)[:, None]) / \
                     np.add(electronic_thermal_conductivity,
                            np.array(lattice_thermal_conductivity)[:, None])

    return zt

def kl(conductivity, seebeck, electronic_thermal_conductivity, zt, temperature):
    """Calculates lattice thermal conductivity.

    Arguments:
        conductivity : array-like
            conductivities.
        seebeck : array-like
            seebeck coefficients.
        electronic_thermal_conductivity : array-like
            electronic thermal conductivities.
        zt : array-like
            zt.
        temperature : array-like
            temperatures.

    Returns:
        np.array
            lattice thermal conductivity.
    """

    pf = power_factor(conductivity, seebeck)
    mid = np.divide(pf * np.array(temperature)[:, None], zt)
    kl = np.subtract(mid, electronic_thermal_conductivity)

    return kl

def power_factor_fromdict(data):
    """Convenience wrapper to calculate power factor from a dictionary.

    Arguments:
        data : dict
            dictionary containing:
                conductivity array-like
                    conductivities.
                seebeck : array-like
                    seebeck coefficients.

    Returns:
        dict
            dictionary with power factors.
    """

    data['power_factor'] = power_factor(data['conductivity'],
                                        data['seebeck'])
    data['meta']['units']['power_factor'] = tp.settings.units()['power_factor']

    return data

def zt_fromdict(data):
    """Convenience wrapper to calculate ZT from a dictionary.

    Arguments:
    data : dict
        dictionary containing:

            conductivity : array-like
                conductivities.
            seebeck : array-like
                seebeck coefficients.
            electronic_thermal_conductivity : array-like
                electronic thermal conductivities.
            lattice_thermal_conductivity : array-like
                lattice thermal conductivities by temperature.
            temperature : array-like
                temperatures in K.

    Returns:
        dict
            dictionary with ZTs.
    """

    data['zt'] = zt(data['conductivity'], data['seebeck'],
                 data['electronic_thermal_conductivity'],
                 data['lattice_thermal_conductivity'],
                 data['temperature'])
    data['meta']['units']['zt'] = tp.settings.units()['zt']

    return data

def kl_fromdict(data):
    """Convenience wrapper to calculate k_latt from a dictionary.

    Arguments:
    data : dict
        dictionary containing:

            conductivity : array-like
                conductivities.
            seebeck : array-like
                seebeck coefficients.
            electronic_thermal_conductivity : array-like
                electronic thermal conductivities.
            zt : array-like
                ZT.
            temperature : array-like
                temperatures in K.

    Returns:
        dict
            dictionary with lattice thermal conductivities.
    """

    q = 'lattice_thermal_conductivity'
    data[q] = kl(data['conductivity'], data['seebeck'],
                 data['electronic_thermal_conductivity'], data['zt'],
                 data['temperature'])
    data['meta']['units'][q] = tp.settings.units()[q]

    return data
