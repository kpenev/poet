"""Create a plot showing the contribution of each tidal term to evolution."""

from os import path

import numpy

from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library
from orbital_evolution.binary import Binary
from orbital_evolution.single_term_interface import SingleTermBody

if __name__ == '__main__':
    e_expansion_fname = path.join(
        path.dirname(path.dirname(path.abspath(__file__))),
        'eccentricity_expansion_coef_O400.sqlite'
    )
    print('Eccentricity expansion: ' + repr(e_expansion_fname))
    assert path.exists(e_expansion_fname)
    orbital_evolution_library.prepare_eccentricity_expansion(
        e_expansion_fname.encode('ascii'),
        1e-4,
        True,
        True
    )
    primary = SingleTermBody(1.0, 1.0)
    primary.set_dissipation(2, 2, 1e-6)
    secondary = SingleTermBody(1.0, 1.0)
    binary = Binary(primary=primary,
                    secondary=secondary,
                    initial_semimajor=1.0,
                    initial_eccentricity=0.3,
                    disk_lock_frequency=1.0,
                    disk_dissipation_age=0.1)
    binary.configure(
        age=1.0,
        semimajor=10.0,
        eccentricity=0.3,
        spin_angmom=(
            (2.0 * numpy.pi / binary.orbital_period(10.0))
            *
            numpy.array([primary.inertia(), secondary.inertia()])
        ),
        inclination=numpy.array([0.0, 0.0]),
        periapsis=numpy.array([0.0]),
        evolution_mode='BINARY'
    )
    print('2,2; 1e-6 Rates: ' + repr(binary.calculate_rates(1.0)))
    primary.set_dissipation(2, 2, 1e-3)
    print('2,2; 1e-3 Rates: ' + repr(binary.calculate_rates(1.0)))
    primary.set_dissipation(1, 0, 1.0)
    print('1,0; 1e0 Rates: ' + repr(binary.calculate_rates(1.0)))
