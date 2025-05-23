"""Settings and defaults.

This module sets the default style sheet, tick locators and axis labels;
as well as providing a means to automatically convert the units
presented and add abbreviations that can be used when loading data.
Custom defaults can be set by saving a copy of tprc.yaml (found in the
main ThermoParser directory) to ``~/.config/tprc.yaml`` and editing
that.
"""

#Functions
#---------
#
#    style:
#        default style sheet.
#    large_style:
#        style sheet for large axes.
#    locator:
#        default tick locators.
#
#
#    to_tp:
#        convert names to tp conventions.
#    to_amset:
#        convert names to amset conventions.
#    to_boltztrap:
#        convert names to boltztrap conventions.
#    to_phono3py:
#        convert names to phono3py conventions.
#
#
#    conversions:
#        default unit conversions.
#    amset_conversions:
#        unit conversions.
#    boltztrap_conversions:
#        unit conversions.
#    boltztrap2_conversions:
#        unit conversions.
#    phonopy_conversions:
#        unit conversions.
#    phono3py_conversions:
#        unit conversions.
#
#
#    units:
#        default units used.
#    dimensions:
#        dimensions of each variable.
#    boltztrap_dimensions:
#        updated dimensions for BoltzTraP.
#
#    labels:
#        default axis labels.
#    inverted_labels:
#        default labels for a dos-style axis.
#    large_labels:
#        default labels for a large axis.
#    long_labels:
#        list of long axis labels.
#    medium_labels:
#        list of slightly abbreviated axis labels.
#    short_labels:
#        list of fully contracted axis labels.
#
#    get_workers:
#        number of workers for parallelisation
#"""

import os
import warnings
import yaml

try:
    filename = '{}/.config/tprc.yaml'.format(os.path.expanduser("~"))
    with open(filename, 'r') as f:
        conf = yaml.safe_load(f)
except yaml.parser.ParserError:
    warnings.warn('Failed to read ~/.config/tprc.yaml')
    conf = None
except FileNotFoundError:
    conf = None

def __dir__():
   """It's a bit of a hack."""

   names = ['locator',
            'to_tp',
            'to_amset',
            'to_boltztrap',
            'to_phono3py',

            'conversions',
            'amset_conversions',
            'boltztrap_conversions',
            'phonopy_conversions',
            'phono3py_conversions',
            'units',
            'dimensions',
            'boltztrap_dimensions',
            'labels',
            'inverted_labels',
            'long_labels',
            'medium_labels',
            'short_labels']

   return names

def style():
    """Get paper-style style sheet."""

    if conf is not None and 'style' in conf and conf['style'] is not None:
        style = conf['style']
    else:
        style = ['tp']
    if isinstance(style, str):
        style = [style]

    return style

def large_style():
    """Get presentation-style style sheet."""

    if conf is not None and 'large_style' in conf and \
       conf['large_style'] is not None:
        style = conf['large_style']
    else:
        style = ['tp-large']
    if isinstance(style, str):
        style = [style]

    return style

def locator():
    """Get default locators."""

    import matplotlib.ticker as ticker


    if conf is not None and 'locator' in conf:
        if 'major' in conf['locator'] and conf['locator']['major'] is not None:
            major = conf['locator']['major']
        else:
            major = 5
        if 'minor' in conf['locator'] and conf['locator']['minor'] is not None:
            minor = conf['locator']['minor']
        else:
            minor = 2
    else:
        major = 5
        minor = 2

    locators = {'major': ticker.MaxNLocator(major),
                'minor': ticker.AutoMinorLocator(minor),
                'log':   ticker.LogLocator(),
                'null':  ticker.NullLocator()}

    return locators

def to_tp():
    """Get dictionary to translate to tp."""

    names = {'ave_pp':               'ph_ph_strength',
             'energies':             'energy',
             'etc':                  'electronic_thermal_conductivity',
             'fermi_levels':         'fermi_level',
             'gv':                   'group_velocity',
             'ir_kpoints':           'kpoint',
             'kpoints':              'kpoint',
             'kappa':                'lattice_thermal_conductivity', # because p3p
             'kappae':               'electronic_thermal_conductivity',
             'kappal':               'lattice_thermal_conductivity',
             'ke':                   'electronic_thermal_conductivity',
             'kl':                   'lattice_thermal_conductivity',
             'ltc':                  'lattice_thermal_conductivity',
             'mfp':                  'mean_free_path',
             'mk':                   'mode_kappa',
             'pf':                   'power_factor',
             'qpoints':              'qpoint',
             'tc':                   'thermal_conductivity',
             'temperatures':         'temperature',
             'scattering_labels':    'stype'}

    if conf is not None and 'to_tp' in conf and conf['to_tp'] is not None:
        names = {**names, **conf['to_tp']}

    return names

def to_amset():
    """Get dictionary to translate to amset."""

    names = {'energy':               'energies',
             'etc':                  'electronic_thermal_conductivity',
             'fermi_level':          'fermi_levels',
             'kappa':                'electronic_thermal_conductivity',
             'kappae':               'electronic_thermal_conductivity',
             'ke':                   'electronic_thermal_conductivity',
             'kpoint':               'ir_kpoints',
             'kpoints':              'ir_kpoints',
             'temperature':          'temperatures',
             'stype':                'scattering_labels'}

    if conf is not None and 'to_amset' in conf and \
       conf['to_amset'] is not None:
        names = {**names, **conf['to_amset']}

    return names

def to_boltztrap():
    """Get dictionary to translate to boltztrap."""

    names = {'etc':                  'electronic_thermal_conductivity',
             'kappa':                'electronic_thermal_conductivity',
             'kappae':               'electronic_thermal_conductivity',
             'ke':                   'electronic_thermal_conductivity'}

    if conf is not None and 'to_boltztrap' in conf and \
       conf['to_boltztrap'] is not None:
        names = {**names, **conf['to_boltztrap']}

    return names

def to_boltztrap2():
    """Get dictionary to translate to boltztrap2."""

    names = {'etc':                  'electronic_thermal_conductivity',
             'kappa':                'electronic_thermal_conductivity',
             'kappae':               'electronic_thermal_conductivity',
             'ke':                   'electronic_thermal_conductivity'}

    if conf is not None and 'to_boltztrap' in conf and \
       conf['to_boltztrap'] is not None:
        names = {**names, **conf['to_boltztrap2']}

    return names

def to_phono3py():
    """Get dictionary to convert to phono3py."""

    names = {'gv':                           'group_velocity',
             'kappal':                       'kappa',
             'kl':                           'kappa',
             'lattice_thermal_conductivity': 'kappa',
             'ltc':                          'kappa',
             'mfp':                          'mean_free_path',
             'mk':                           'mode_kappa',
             'ph_ph_strength':               'ave_pp',
             'qpoint':                       'qpoint',
             'temperatures':                 'temperature'}

    if conf is not None and 'to_phono3py' in conf and \
       conf['to_phono3py'] is not None:
        names = {**names, **conf['to_phono3py']}

    return names

def conversions():
    """Get dictionary of custom unit conversions from tprc.yaml."""

    conversions = {}

    if conf is not None and 'conversions' in conf and \
       conf['conversions'] is not None:
        conversions = conf['conversions']

    return conversions

def amset_conversions():
    """Get dictionary of unit conversions for amset to tp."""

    conversions = {}

    return conversions

def boltztrap_conversions():
    """Get dictionary of unit conversions for boltztrap to tp."""

    conversions = {}

    return conversions

def boltztrap2_conversions():
    """Get dictionary of unit conversions for boltztrap to tp."""

    conversions = {}

    return conversions

def phonopy_conversions():
    """Get dictionary of unit conversions for phonopy to tp."""

    conversions = {}

    return conversions

def phono3py_conversions():
    """Get dictionary of unit conversions for phono3py to tp."""

    conversions = {'group_velocity': 1e2,         # THz AA to m s-1
                   'gv_by_gv':       1e4,         # THz2 AA2 to m2 s-2
                   'heat_capacity':  1.60218e-19, # eV K-1 to J K-1
                   'lifetime':       1e-12,       # THz-1 to s
                   'mean_free_path': 1e-10}       # AA to m

    return conversions

def units(use_tprc=True):
    """Get dictionary of units of quantities used in tp.

    Arguments
    ---------

        use_tprc : bool, optional
            read custom units from tprc.yaml if present. Default: True.
    """

    units = {'average_eff_mass':                'm_e',
             'chemical_potential':              'eV',
             'conductivity':                    'S m-1',
             'doping':                          'cm-3',
             'efermi':                          'eV',
             'effective_mass':                  'm_e',
             'electronic_heat_capacity':        'J mol-1 K-1',
             'electronic_thermal_conductivity': 'W m-1 K-1',
             'energy':                          'eV',
             'fd_weights':                      '',
             'fermi_level':                     'eV',
             'frequency':                       'THz',
             'gamma':                           'THz',
             'group_velocity':                  'm s-1',
             'gruneisen':                       '',
             'gv_by_gv':                        'm2 s-2',
             'hall':                            'm3 C-1',
             'hall_carrier_concentration':      'cm-3',
             'heat_capacity':                   'J K-1',
             'ibz_weights':                     '',
             'kpoint':                          '',
             'lattice_thermal_conductivity':    'W m-1 K-1',
             'lifetime':                        's',
             'mean_free_path':                  'm',
             'mobility':                        'cm^2 V-1 s-1',
             'mode_kappa':                      'W m-1 K-1',
             'mu_bounds':                       'eV',
             'normalised_weights':              '',
             'occupation':                      'phonons',
             'ph_ph_strength':                  'eV2',
             'power_factor':                    'W m-1 K-2',
             'qpoint':                          '',
             'scattering_rates':                's-1',
             'seebeck':                         'muV K-1',
             'seebeck_effective_mass':          'm_e',
             'temperature':                     'K',
             'thermal_conductivity':            'W m-1 K-1',
             'total_weights':                   '',
             'velocities':                      'm s-1',
             'weighted_mfp':                    'm',
             'weighted_rates':                  's-1',
             'zt':                              ''}

    if use_tprc and conf is not None and 'units' in conf and \
       conf['units'] is not None:
        units = {**units, **conf['units']}

    return units

def dimensions():
    """Get dictionary of dimensions of quantities used in tp."""

    dims = {'average_eff_mass':                ['temperature', 'doping', 3, 3],
            'chemical_potential':              [],
            'conductivity':                    ['temperature', 'doping', 3, 3],
            'doping':                          ['doping'],
            'efermi':                          [],
            'effective_mass':                  ['temperature', 'doping', 3, 3],
            'electronic_heat_capacity':        ['temperature', 'doping'],
            'electronic_thermal_conductivity': ['temperature', 'doping', 3, 3],
            'energy':                          ['band', 'kpoint'],
            'fd_weights':                      ['temperature', 'doping', 'band', 'kpoint'],
            'fermi_level':                     ['temperature', 'doping'],
            'frequency':                       ['qpoint', 'band'],
            'gamma':                           ['temperature', 'qpoint', 'band'],
            'group_velocity':                  ['qpoint', 'band', 3],
            'gruneisen':                       ['qpoint', 'band'],
            'gv_by_gv':                        ['qpoint', 'band', 6],
            'hall':                            ['temperature', 'doping', 3, 3, 3],
            'hall_carrier_concentration':      [],
            'heat_capacity':                   ['temperature', 'qpoint', 'band'],
            'ibz_weights':                     ['kpoint'],
            'kpoint':                          ['kpoint', 3],
            'lattice_thermal_conductivity':    ['temperature', 6],
            'lifetime':                        ['temperature', 'qpoint', 'band'],
            'mean_free_path':                  ['temperature', 'qpoint', 'band', 3],
            'mesh':                            [3],
            'mobility':                        ['stype', 'temperature', 'doping', 3, 3],
            'mode_kappa':                      ['temperature', 'qpoint', 'band', 6],
            'mu_bounds':                       [],
            'normalised_weights':              ['temperature', 'doping', 'band', 'kpoint'],
            'occupation':                      ['temperature', 'qpoint', 'band'],
            'ph_ph_strength':                  ['qpoint', 'band'],
            'power_factor':                    ['temperature', 'doping', 3, 3],
            'qpoint':                          ['qpoint', 3],
            'scattering_rates':
                ['stype', 'doping', 'temperature', 'band', 'kpoint'],
            'seebeck':                         ['temperature', 'doping', 3, 3],
            'seebeck_effective_mass':          [],
            'temperature':                     ['temperature'],
            'thermal_conductivity':            ['temperature', 'doping', 3, 3],
            'total_weights':                   ['temperature', 'doping'],
            'velocities':                      ['band', 'kpoint', 3],
            'weighted_mfp':                    ['stype', 'temperature', 'doping', 3],
            'weighted_rates':                  ['stype', 'temperature', 'doping'],
            'zt':                              ['temperature', 'doping', 3, 3]}

    return dims

def boltztrap_dimensions():
    """Get dictionary of dimensions of quantities updated for BoltzTraP."""

    dims = dimensions()
    dims['average_eff_mass'] =                ['dtype', 'temperature', 'doping', 3, 3]
    dims['conductivity'] =                    ['dtype', 'temperature', 'doping', 3, 3]
    dims['electronic_thermal_conductivity'] = ['dtype', 'temperature', 'doping', 3, 3]
    dims['fermi_level'] =                     ['dtype', 'temperature', 'doping']
    dims['mobility'] =                        ['dtype', 'temperature', 'doping', 3, 3]
    dims['power_factor'] =                    ['dtype', 'temperature', 'doping', 3, 3]
    dims['seebeck'] =                         ['dtype', 'temperature', 'doping', 3, 3]

    return dims

def labels():
    """Get the default labels for small axes in tp."""

    labels = {'short':  short_labels,
              'medium': medium_labels,
              'long':   long_labels}

    if conf is not None and 'labels' in conf and conf['labels'] is not None:
        length = conf['labels']
    else:
        length = 'long'

    return labels[length]()

def inverted_labels():
    """Get the default labels for inverted axes in tp."""

    labels = {'short':  short_labels,
              'medium': medium_labels,
              'long':   long_labels}

    if conf is not None and 'inverted_labels' in conf and \
       conf['inverted_labels'] is not None:
        length = conf['inverted_labels']
    else:
        length = 'short'

    return labels[length]()

def large_labels():
    """Get the default labels for large axes."""

    labels = {'short':  short_labels,
              'medium': medium_labels,
              'long':   long_labels}

    if conf is not None and 'large_labels' in conf and \
       conf['large_labels'] is not None:
        length = conf['large_labels']
    else:
        length = 'medium'

    return labels[length]()

def long_labels():
    """Get a dictionary of long-form axis labels."""

    labels = {'chemical_potential':
                  'Chemical Potential (eV)',
              'complexity_factor':
                  'Complexity Factor',
              'conductivity':
                  r'Conductivity (S m$\mathregular{^{-1}}$)',
              'cumulative_kappa':
                  r'Cumulative Lattice Thermal Conductivity (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'cumulative_percent':
                  'Cumulative Lattice Thermal Conductivity (%)',
              'doping':
                  r'Carrier Concentration (cm$\mathregular{^{-3}}$)',
              'dos':
                  'Density of States',
              'efermi':
                  'Fermi Energy (eV)',
              'effective_mass':
                  r'Effective Mass (m$\mathregular{_e}$)',
              'energy':
                  'Energy (eV)',
              'electronic_heat_capacity':
                  r'Electronic Specific Heat Capacity (J mol$\mathregular{^{-1}\ K^-1}$)',
              'electronic_thermal_conductivity':
                  r'Electronic Thermal Conductivity (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'fermi_level':
                  'Fermi Level (eV)',
              'frequency':
                  'Frequency (THz)',
              'gamma':
                  'Imaginary Self Energy (THz)',
              'group_velocity':
                  r'Group Velocity (m s$\mathregular{^{-1}}$)',
              'gruneisen':
                  'Gruneisen Parameter',
              'gv_by_gv':
                  r'Group Velocity Outer Product (m$\mathregular{^2\ s^{-2}}$)',
              'hall':
                  r'Hall Coefficient (m$\mathregular{^3\ C^{-1}}$)',
              'heat_capacity':
                  r'Heat Capacity (J K$\mathregular{^{-1}}$)',
              'lattice_thermal_conductivity':
                  r'Lattice Thermal Conductivity (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'lifetime':
                  'Lifetime (s)',
              'mean_free_path':
                  'Mean Free Path (m)',
              'mobility':
                  r'Mobility (cm$\mathregular{^2\ V^{-1}\ s^{-1}}$)',
              'mode_kappa':
                  r'Lattice Thermal Conductivity (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'ph_ph_strength':
                  r'Average Phonon-Phonon Interaction Strengths (eV$\mathregular{^2}$)',
              'power_factor':
                  r'Power Factor (W m$\mathregular{^{-1}\ K^{-2}}$)',
              'occupation':
                  'Occupation',
              'scattering_rates':
                  r'Scattering Rates (s$\mathregular{^{-1}}$)',
              'seebeck':
                  r'Seebeck Coefficient ($\mathregular{\mu V\ K^{-1}}$)',
              'temperature':
                  'Temperature (K)',
              'thermal_conductivity':
                  r'Thermal Conductivity (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'velocities':
                  r'Velocity (m s$\mathregular{^{-1}})',
              'wavevector':
                  'Wavevector',
              'weighted_mfp':
                  'Mean Free Path (m)',
              'weighted_rates':
                  r'Scattering Rates (s$\mathregular{^{-1}}$)',
              'zt':
                  'ZT'}

    if conf is not None and 'long_labels' in conf and \
       conf['long_labels'] is not None:
        labels = {**labels, **conf['long_labels']}

    return labels

def medium_labels():
    """Get a dictionary of medium-length axis labels."""

    labels = {'chemical_potential':
                  'Chemical Potential (eV)',
              'complexity_factor':
                  'Complexity Factor',
              'conductivity':
                  r'Conductivity (S m$\mathregular{^{-1}}$)',
              'cumulative_kappa':
                  r'Cum. LTC (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'cumulative_percent':
                  'Cum. LTC (%)',
              'doping':
                  r'Carrier Concentration (cm$\mathregular{^{-3}}$)',
              'dos':
                  'Density of States',
              'efermi':
                  'Fermi Energy (eV)',
              'effective_mass':
                  r'Effective Mass (m$\mathregular{_e}$)',
              'energy':
                  'Energy (eV)',
              'electronic_heat_capacity':
                  r'Elec. Spec. Heat Cap. (J mol$\mathregular{^{-1}\ K^-1}$)',
              'electronic_thermal_conductivity':
                  r'Elec. Therm. Cond. (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'fermi_level':
                  'Fermi Level (eV)',
              'frequency':
                  'Frequency (THz)',
              'gamma':
                  'Imaginary Self Energy (THz)',
              'group_velocity':
                  r'Group Velocity (m s$\mathregular{^{-1}}$)',
              'gruneisen':
                  'Gruneisen Parameter',
              'gv_by_gv':
                  r'Group Vel. Outer Prod. (m$\mathregular{^2\ s^{-2}}$)',
              'hall':
                  r'Hall Coefficient (m$\mathregular{^3\ C^{-1}}$)',
              'heat_capacity':
                  r'Heat Capacity (J K$\mathregular{^{-1}}$)',
              'lattice_thermal_conductivity':
                  r'Lat. Therm. Cond. (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'lifetime':
                  'Lifetime (s)',
              'mean_free_path':
                  'Mean Free Path (m)',
              'mobility':
                  r'Mobility (cm$\mathregular{^2\ V^{-1}\ s^{-1}}$)',
              'mode_kappa':
                  r'Lat. Therm. Cond. (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'ph_ph_strength':
                  r'Avg. Ph-Ph Strengths (eV$\mathregular{^2}$)',
              'power_factor':
                  r'Power Factor (W m$\mathregular{^{-1}\ K^{-2}}$)',
              'occupation':
                  'Occupation',
              'scattering_rates':
                  r'Scattering Rates (s$\mathregular{^{-1}}$)',
              'seebeck':
                  r'Seebeck Coefficient ($\mathregular{\mu V\ K^{-1}}$)',
              'temperature':
                  'Temperature (K)',
              'thermal_conductivity':
                  r'Thermal Cond. (W m$\mathregular{^{-1}\ K^{-1}}$)',
              'velocities':
                  r'Velocity (m s$\mathregular{^{-1}})',
              'wavevector':
                  'Wavevector',
              'weighted_mfp':
                  'Mean Free Path (m)',
              'weighted_rates':
                  r'Scattering Rates (s$\mathregular{^{-1}}$)',
              'zt':
                  'ZT'}

    if conf is not None and 'medium_labels' in conf and \
       conf['medium_labels'] is not None:
        labels = {**labels, **conf['medium_labels']}

    return labels

def short_labels():
    """Get dictionary of short-form axis labels."""

    labels = {'chemical_potential':
                  r'$\mathregular{\mu}$ (eV)',
              'complexity_factor':
                  r'$\mathregular{N_v*K*}$',
              'conductivity':
                  r'$\mathregular{\sigma\ (S\ m^{-1})}$',
              'cumulative_kappa':
                  r'$\mathregular{\kappa_l\ (W\ m^{-1}\ K^{-1})}$',
              'cumulative_percent':
                  r'$\mathregular{\kappa_l}$ (%)',
              'doping':
                  r'n (cm$\mathregular{^{-3}}$)',
              'dos':
                  'DoS',
              'efermi':
                  r'E$\mathregular{_{F}}$ (eV)',
              'effective_mass':
                  r'$\mathregular{m*\ (m_e})$',
              'energy':
                  'E (eV)',
              'electronic_heat_capacity':
                  r'e- c_V (J mol$\mathregular{^{-1}\ K^-1}$)',
              'electronic_thermal_conductivity':
                  r'$\mathregular{\kappa_e\ (W\ m^{-1}\ K^{-1})}$',
              'fermi_level':
                  r'E$\mathregular{_{F}}$ (eV)',
              'frequency':
                  r'$\mathregular{\\nu}$ (THz)',
              'gamma':
                  r'$\mathregular{\Gamma}$ (THz)',
              'group_velocity':
                  r'$\mathregular{g_v\ (m\ s^{-1})}$',
              'gruneisen':
                  r'$\gamma$',
              'gv_by_gv':
                  r'$\mathregular{g_v \otimes g_v\ (m^2\ s^{-2})}$',
              'hall':
                  r'$\mathregular{R_H\ (m^3\ C^{-1}}$)',
              'heat_capacity':
                  r'C_V (J K$\mathregular{^{-1}}$)',
              'lattice_thermal_conductivity':
                  r'$\mathregular{\kappa_l\ (W\ m^{-1}\ K^{-1})}$',
              'lifetime':
                  r'$\mathregular{\\tau}$ (s)',
              'mean_free_path':
                  r'$\mathregular{\Lambda}$ (m)',
              'mobility':
                  r'$\mathregular{\mu\ (cm^2\ V^{-1}\ s^{-1})}$',
              'mode_kappa':
                  r'$\mathregular{\kappa_l\ (W\ m^{-1}\ K^{-1})}$',
              'ph_ph_strength':
                  r'$\mathregular{P_\lambda\ (eV^2)}$',
              'power_factor':
                  r'PF (W m$\mathregular{^{-1}\ K^{-2}}$)',
              'occupation':
                  'Occupation',
              'scattering_rates':
                  r'$\mathregular{\\tau^{-1}\ (s^{-1})}$',
              'seebeck':
                  r'$\mathregular{\\alpha\ (\mu V\ K^{-1})}$',
              'temperature':
                  'T (K)',
              'thermal_conductivity':
                  r'$\mathregular{\kappa\ (W\ m^{-1}\ K^{-1}}$)',
              'velocities':
                  r'v (m s$\mathregular{^{-1}})',
              'wavevector':
                  'q',
              'weighted_mfp':
                  r'$\mathregular{\Lambda}$ (m)',
              'weighted_rates':
                  r'$\mathregular{\\tau^{-1}\ (s^{-1})}$',
              'zt':
                  'ZT'}

    if conf is not None and 'short_labels' in conf and \
       conf['short_labels'] is not None:
        labels = {**labels, **conf['short_labels']}

    return labels

def get_workers():
    """Get default number of workers for parallelisation.
    
    Defaults to number of cores.
    """

    workers = os.cpu_count()

    if conf is not None and 'workers' in conf and \
       conf['workers'] is not None:
        workers = int(conf['workers'])

    return workers
