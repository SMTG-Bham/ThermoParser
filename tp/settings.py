"""Settings for appearance, naming, units and labels.

Also adds abbreviations. In future may have to split up conversions.

Functions:
    style:
        default style sheet.
    locator:
        default tick locators.
    to_tp:
        convert names to tp conventions.
    to_amset:
        convert names to amset conventions.
    to_amset:
        convert names to boltztrap conventions.
    to_phono3py:
        convert names to phono3py conventions.
    conversions:
        unit conversions (see tp.data.load).
    units:
        default units used.
    labels:
        default axis labels.
"""

def __dir__():
    names = ['to_tp',
             'to_amset',
             'to_boltztrap',
             'to_phono3py',

             'conversions',
             'units',
             'labels']

    return names

def style():
    """Get default style sheet."""
    return 'tp'

def locator():
    """Get default locators"""

    import matplotlib.ticker as ticker

    locators = {'major': ticker.MaxNLocator(5),
                'minor': ticker.AutoMinorLocator(2),
                'log':   ticker.LogLocator()}

    return locators

def to_tp():
    """Get dictionary to translate to tp."""

    names = {'energies':     'energy',
             'fermi_levels': 'fermi_level',
             'gv':           'group_velocity',
             'kappa':        'lattice_thermal_conductivity', # because p3p
             'kappae':       'electronic_thermal_conductivity',
             'kappal':       'lattice_thermal_conductivity',
             'ke':           'electronic_thermal_conductivity',
             'kl':           'lattice_thermal_conductivity',
             'mfp':          'mean_free_path',
             'mk':           'mode_kappa',
             'pf':           'power_factor',
             'temperatures': 'temperature'}

    return names

def to_amset():
    """Get dictionary to translate to amset."""

    names = {'energy':      'energies',
             'fermi_level': 'fermi_levels',
             'kappa':       'electronic_thermal_conductivity',
             'kappae':      'electronic_thermal_conductivity',
             'ke':          'electronic_thermal_conductivity',
             'temperature': 'temperatures'}

    return names

def to_boltztrap():
    """Get dictionary to translate to boltztrap."""

    names = {'electronic_thermal_conductivity': 'thermal_conductivity',
             'kappa':                           'thermal_conductivity',
             'kappae':                          'thermal_conductivity',
             'ke':                              'thermal_conductivity'}

    return names

def to_phono3py():
    """Get dictionary to convert to phono3py variable names."""

    names = {'gv':                           'group_velocity',
             'kappal':                       'kappa',
             'kl':                           'kappa',
             'lattice_thermal_conductivity': 'kappa',
             'mfp':                          'mean_free_path',
             'mk':                           'mode_kappa',
             'temperatures':                 'temperature'}

    return names

def amset_conversions():
    """Get dictionary of unit conversions for phonopy."""

    conversions = {}

    return conversions

def boltztrap_conversions():
    """Get dictionary of unit conversions for phonopy."""

    conversions = {}

    return conversions

def phonopy_conversions():
    """Get dictionary of unit conversions for phonopy."""

    conversions = {'power_factor': 1e-6} # muW m-1 K-2 to W m-1 K-2

    return conversions

def phono3py_conversions():
    """Get dictionary of unit conversions for phono3py."""

    conversions = {'group_velocity': 1e2,         # THz AA to m s-1
                   'gv_by_gv':       1e4,         # THz2 AA2 to m2 s-2
                   'heat_capacity':  1.60218e-19, # eV K-1 to J K-1
                   'lifetime':       1e-12,       # THz-1 to s
                   'mean_free_path': 1e-10}       # AA to m

    return conversions

def units():
    """Get dictionary of units of quantities used in tp"""

    units = {'average_eff_mass':                'm_e',
             'chemical_potential':              'eV',
             'conductivity':                    'S m-1',
             'doping':                          'cm-3',
             'efermi':                          'eV',
             'effective_mass':                  'm_e',
             'electronic_thermal_conductivity': 'W m-1 K-1',
             'energy':                          'eV',
             'fermi_level':                     'eV',
             'frequency':                       'THz',
             'gamma':                           'THz',
             'group_velocity':                  'm s-1',
             'gv_by_gv':                        'm2 s-2',
             'hall_carrier_concentration':      'cm-3',
             'heat_capacity':                   'J K-1',
             'lattice_thermal_conductivity':    'W m-1 K-1',
             'lifetime':                        's',
             'mean_free_path':                  'm',
             'mobility':                        'cm^2 V-1 s-1',
             'mode_kappa':                      'W m-1 K-1',
             'mu_bounds':                       'eV',
             'occupation':                      'phonons',
             'power_factor':                    'W m-1 K-2',
             'scattering_rates':                's-1',
             'seebeck':                         'muV K-1',
             'seebeck_effective_mass':          'm_e',
             'temperature':                     'K',
             'zt':                              ''}

    return units

def labels():
    """Get dictionary of axis labels used in tp"""

    labels = {'chemical_potential':
                  'Chemical Potential (eV)',
              'complexity_factor':
                  'Complexity Factor',
              'conductivity':
                  'Conductivity (S m$\mathregular{^{-1}}$)',
              'cumulative_kappa':
                  'Cumulative Lattice Thermal Conductivity (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'doping':
                  'Carrier Concentration (cm$\mathregular{^{-1}}$)',
              'dos':
                  '',#'Arbitrary Units',#
              'efermi':
                  'Fermi Energy (eV)',
              'effective_mass':
                  'Effective Mass (m$\mathregular{_e}$)',
              'energy':
                  'Energy (eV)',
              'electronic_thermal_conductivity':
                  'Electronic Thermal Conductivity (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'fermi_level':
                  'Fermi Level (eV)',
              'frequency':
                  'Frequency (THz)',
              'gamma':
                  'Imaginary Self Energy (THz)',
              'group_velocity':
                  'Group Velocity (m s$\mathregular{^{-1}}$)',
              'gv_by_gv':
                  'Group Velocity Outer Product (m$\mathregular{^2\ s^{-2}}$)',
              'heat_capacity':
                  'Heat Capacity (J K$\mathregular{^{-1}}$)',
              'lattice_thermal_conductivity':
                  'Lattice Thermal Conductivity (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'lifetime':
                  'Lifetime (s)',
              'mean_free_path':
                  'Mean Free Path (m)',
              'mobility':
                  'Mobility (cm$\mathregular{^2\ V^{-1}\ s^{-1}}$)',
              'mode_kappa':
                  'Lattice Thermal Conductivity (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'power_factor':
                  'Power Factor (W m$\mathregular{^{-1}\ K^{-2}}$)',
              'occupation':
                  'Occupation',
              'scattering_rates':
                  'Scattering Rates (s$\mathregular{^{-1}}$)',
              'seebeck':
                  'Seebeck Coefficient ($\mathregular{\mu V\ K^{-1}}$)',
              'temperature':
                  'Temperature (K)',
              'wavevector':
                  'Wavevector',
              'zt':
                  'ZT'}

    return labels
