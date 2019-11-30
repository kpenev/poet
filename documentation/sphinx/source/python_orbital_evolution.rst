***********************************
The Python Orbital Evolution Module
***********************************

A module for calculating orbital evolutions of two body systems in cilcular
orbits. Both bodies are split into a number of zones, some or all of which are
dissipative. The evolution tracks each zone's angular momentum vector allowed to
point in an arbitrary direction. All neighboring zones exchange angular momentum
and the outermost zone possibly loses angular momentum due to magnetically
launched wind. All frequency components of the dissipation are treated
separately with their own tidal frequency. The dissipation efficiency is a
broken powerlaw of tidal and stellar spin frequency.

A binray system in poet consists of a primary and a secondary object. With the
primary object forming first and evolving in isolation until a user specified
time when the secondary object appears in the system, and binary evolution
commences. The secondary must be configured before running the evolution to the
state in which it will first appear. Two types of objects can be part of a
POET binary system: stars and planets.

Planets in POET are quite simple. They consist of a single zone of fixed mass,
radius and moment of inertia, and do not lose angular momentum to wind.
Dissipation in the (only zone of the) planet can be specified the same way it is
specified for any other zone (see below).

With the default stellar evolution included with POET, stars must satisfy:
:math:`0.4 M_\odot \leq M_\star \leq 1.2 M_\odot`, and :math:`-1 < [Fe/H] <
0.5`. The have two zones: a radiative core and a convective envelope, with the
boundary betweeen these zones both in mass and in radius evolving with age. The
dissipation of each zone can be specified independently, though theoretically
the expectation is that the dissipation in the convective zone will have a much
much bigger effect than any likely dissipation in the core. In addition, the
zones exchange angular momentum, and the radiative envelope loses angular
momentum to stelar wind.

In order to calculate the evolution of a planet-star system the basic steps are:

  * Create a stellar evolution interpolator

  * Create, configure and initialize a star

  * Create a planet

  * Combine the two objects in a binary

  * Calculate the evolution and query the results

Examples on how to accomplish each step follow. Note that if a zone of an object
is not dissipative it is better to avoid specifying the disispation rather than
specifying a zero phase-lag, since in the latter case tidal torque and power are
still computed though they are guaranteed to always be zero, which just wastes
CPU time:

Create a Stellar Evolution Interpolator
=======================================

More details about stellar evolution interpolators are given :doc:`here
<python_stellar_evolution>`, but the basic steps are::

    from stellar_evolution.manager import StellarEvolutionManager
    import numpy

    #Create a manager
    manager = StellarEvolutionManager(<path to serialized interpolators>)

    #Find and use the default interpolator
    interpolator = manager.get_interpolator_by_name('default')

Create, Configure and Initialize a Star
=======================================

Given an interpolator a star with a given mass, [Fe/H], spin-down model
parameters and dissipation parameters is created and prepared for evolution as
follows::

    from orbital_evolution.transformations import phase_lag
    from orbital_evolution.star_interface import EvolvingStar

    #Create the star with
    #   mass of 1 solar mass
    #   [Fe/H] = 0
    #   Spin-down parameters: K_w = 0.17 in solar units times rad/day per Gyr
    #                         \omega_{sat} = 2.45 rad/day
    #   core-envelope coupling timescale of 5Myrs
    star = EvolvingStar(mass = 1.0,
                        metallicity = 0.0,
                        wind_strength = 0.17,
                        wind_saturation_frequency = 2.45,
                        diff_rot_coupling_timescale = 5.0e-3,
                        interpolator = interpolator)

    #Prepare the stellar evolution interpolation to start following the
    #evolution of the star as soon as the core stars to form.
    star.select_interpolation_region(star.core_formation_age())

    #Define the dissipation of the stellar envelope. In this case, a constant
    #phase lag corresponding to \f$Q'_\star = 10^6\f$.
    star.set_dissipation(zone_index = 0,
                         tidal_frequency_breaks = None,
                         spin_frequency_breaks = None,
                         tidal_frequency_powers = numpy.array([0.0]),
                         spin_frequency_powers = numpy.array([0.0]),
                         reference_phase_lag = phase_lag(6.0))

    #The dissipation in the core of the star is not specified which results in a
    #non-disspative core

Create a Planet
===============

Creating planets is very simple. One only needs to specify the mass and radius
**in soler units**, we do not specify dissipation, resulting in a
non-dissipative planet::

    from orbital_evolution.planet_interface import LockedPlanet
    from astropy import units, constants

    #Create a planet with a Jovian mass and radius
    planet = LockedPlanet(
        mass = constants.M_jup.to_value('M_sun'),
        radius = constants.R_jup.to_value('R_sun')
    )

Combine the Two Objects in a Binary
===================================

A star and a planet can be combined in a binary which is then prepared for
evolution calculations as follows::

    from orbital_evolution.binary import Binary
    from math import pi
    import numpy

    #Create a binary with the primary object being the star and the seconday
    #being an initially non-existent the planet, which forms when the
    #protoplanetary disk dissipates (in this case 4Myrs) with an initial orbital
    #period of 3 days in a circular and aligned orbit.
    #Prior to the planet formith the stellar envelope spin frequency is held
    #fixed (presumably locked to the inner edge of the disk) at a period of 7
    #days.
    disk_dissipation_age = 4e-3
    porb_initial = 3.0
    binary = Binary(primary = star,
                    secondary = planet,
                    initial_orbital_period = porb_initial,
                    initial_eccentricity = 0.0,
                    initial_inclination = 0.0,
                    disk_lock_frequency = 2.0 * pi / 7.0,
                    disk_dissipation_age = disk_dissipation_age,
                    secondary_formation_age = 4e-3)

    #Specify the initial conditions at which the binary starts. In thisa cese
    #trivial since the planet does not exist at the start of the evolution
    #calculations.
    binary.configure(age = star.core_formation_age(),
                     semimajor = float('nan'),
                     eccentricity = float('nan'),
                     spin_angmom = numpy.array([0.0]),
                     inclination = None,
                     periapsis = None,
                     evolution_mode = 'LOCKED_SURFACE_SPIN')

    #Configure the planet in the state in which it will form.
    planet.configure(age = disk_dissipation_age,
                     companion_mass = star.mass,
                     semimajor = binary.semimajor(porb_initial),
                     eccentricity = 0.0,
                     spin_angmom = numpy.array([0.0]),
                     inclination = None,
                     periapsis = None,
                     locked_surface = False,
                     zero_outer_inclination = True,
                     zero_outer_periapsis = True)

    #Configuring the binary (the first command above) will trigger configuring
    #the star, and before evolution starts the star needs to set its wind
    #saturation flag per its current state.
    star.detect_stellar_wind_saturation()

Calculate the evolution and query the results
=============================================

Once we have a binary system, the following calculates the evolution and obtains
the specified quantities at each timestep::

    from orbital_evolution.evolve_interface import library as\
    orbital_evolution_library

    #Read-in the pre-calculated coefficients giving the expansion of the tidal
    #potential as a series in eccentricity and frequency.
    #The filename MUST BE SUPPLIED AS BYTES (not string)!!!
    orbital_evolution_library.read_eccentricity_expansion_coefficients(
        <filename with eccentricity expansion coefficients>
    )

    #Calculate the evolution to a final age of 10Gyrs
    #With a maximum step size of 1Myr
    #With the ODE solver precision set ot 10^-6
    #An empty list of ages which the evolution must precisely visit
    binary.evolve(10.0, 0.001, 1e-6, None)

    #Select the quantities to get at each evolution timestep.
    #Further requests for quantities can be made later without re-calculating
    #the evolution.
    evolution_quantities = ['age',
                            'semimajor',
                            'envelope_angmom',
                            'core_angmom',
                            'wind_saturation']
    
    #Get the specified quantities as a numpy record array with record names
    #given by the quantity names above.
    evolution = binary.get_evolution(evolution_quantities)
