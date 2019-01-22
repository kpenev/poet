"""Define class for finding initial conditions to match current properties."""

from scipy.optimize import brentq
import scipy

from orbital_evolution.binary import Binary
from orbital_evolution.star_interface import EvolvingStar
from basic_utils import Structure

#TODO: organize attributes into dictionaries or structures.
#pylint: disable=too-many-instance-attributes
class InitialConditionSolver:
    """Find initial conditions which reproduce a given system now."""

    def _try_initial_conditions(self, initial_orbital_period, disk_period):
        """
        Get present orbital and stellar spin periods for initial conditions.

        Args:
            initial_orbital_period:    The initial orbital period to calculate
                the deviation for.

            disk_period:    The disk locking period to calculate the deviation
                for.

        Returns:
            (float, float):
                The present day orbital period of the system resulting when
                the evolution is started with the input periods.

                The present day surface spin of the star resulting when the
                evolution is started with the input periods.
        """

        print('Trying P0 = %s, Pdisk = %s'
              %
              (repr(initial_orbital_period), repr(disk_period)))
        if self.binary is not None:
            self.binary.delete()
        self.binary = Binary(
            primary=self.primary,
            secondary=self.secondary,
            initial_orbital_period=initial_orbital_period,
            initial_eccentricity=self.initial_eccentricity,
            initial_inclination=0.0,
            disk_lock_frequency=2.0 * scipy.pi / disk_period,
            disk_dissipation_age=self.target.disk_dissipation_age,
            secondary_formation_age=self.target.planet_formation_age
        )
        print('Constructed binary')

        self.primary.select_interpolation_region(
            self.primary.core_formation_age()
        )

        self.binary.configure(age=self.primary.core_formation_age(),
                              semimajor=float('nan'),
                              eccentricity=float('nan'),
                              spin_angmom=scipy.array([0.0]),
                              inclination=None,
                              periapsis=None,
                              evolution_mode='LOCKED_SURFACE_SPIN')
        print('Configured binary')

        self.primary.detect_stellar_wind_saturation()
        if isinstance(self.secondary, EvolvingStar):
            secondary_inclination_periapsis = dict(
                inclination=scipy.array([0.0]),
                periapsis=scipy.array([0.0])
            )
        else:
            secondary_inclination_periapsis = dict(
                inclination=None,
                periapsis=None
            )
        self.secondary.configure(
            age=self.target.planet_formation_age,
            companion_mass=self.binary.primary.mass,
            semimajor=self.binary.semimajor(initial_orbital_period),
            eccentricity=self.initial_eccentricity,
            spin_angmom=scipy.array([0.0, 0.0]),
            locked_surface=False,
            zero_outer_inclination=True,
            zero_outer_periapsis=True,
            **secondary_inclination_periapsis
        )
        print('Configured secondary')
        if isinstance(self.secondary, EvolvingStar):
            self.secondary.detect_stellar_wind_saturation()
        self.binary.evolve(
            self.target.age,
            self.evolution_max_time_step,
            self.evolution_precision,
            None,
            True
        )
        final_state = self.binary.final_state()
        #False positives
        #pylint: disable=no-member
        assert final_state.age == self.target.age
        orbital_period = self.binary.orbital_period(final_state.semimajor)
        stellar_spin_period = (
            2.0 * scipy.pi
            *
            self.binary.primary.envelope_inertia(final_state.age)
            /
            final_state.primary_envelope_angmom
        )
        #pylint: enable=no-member
        print('Got Porb = %s, P* = %s'
              %
              (repr(orbital_period), repr(stellar_spin_period)))

        if scipy.isnan(orbital_period):
            orbital_period = 0.0
        return orbital_period, stellar_spin_period

    def _find_porb_range(self, guess_porb_initial, disk_period):
        """
        Find initial orbital period range where final porb error flips sign.

        Args:
            guess_porb_initial:    An initial guess for where the sign change
                occurs.

            disk_period:    The disk locking period to assume during the search.

        Returns:
            (float, float):
                A pair of initial orbital periods for which the sign of the
                final orbital period error changes.
        """

        porb_min, porb_max = scipy.nan, scipy.nan
        porb_initial = guess_porb_initial
        porb = self._try_initial_conditions(porb_initial, disk_period)[0]
        porb_error = porb - self.target.Porb
        guess_porb_error = porb_error
        step = 2.0 if guess_porb_error < 0 else 0.5

        while porb_error * guess_porb_error > 0 and porb_initial < 100.0:
            if porb_error < 0:
                porb_min = porb_initial
            else:
                porb_max = porb_initial
            porb_initial *= step
            porb = self._try_initial_conditions(porb_initial, disk_period)[0]
            if not scipy.isnan(porb):
                porb_error = porb - self.target.Porb

        if scipy.isnan(porb_error):
            return scipy.nan, scipy.nan

        if porb_error < 0:
            porb_min = porb_initial
        else:
            porb_max = porb_initial
            if porb_error == 0:
                porb_min = porb_initial

        print('For Pdisk = %s, orbital period range: %s < Porb < %s'
              %
              (repr(disk_period), repr(porb_min), repr(porb_max)))
        return porb_min, porb_max

    def __init__(self,
                 *,
                 planet_formation_age=None,
                 disk_dissipation_age=None,
                 evolution_max_time_step=1.0,
                 evolution_precision=1e-6,
                 orbital_period_tolerance=1e-6,
                 spin_tolerance=1e-6,
                 initial_eccentricity=0.0):
        """
        Initialize the object.

        Args:
            planet_formation_age:    If not None, the planet is assumed to form
                at the given age (in Gyr). Otherwise, the starting age must be
                specified each time this object is called.

            disk_dissipation_age:    The age at which the disk dissipates in
                Gyrs.

            evolution_max_time_step:    The maximum timestep the evolution is
                allowed to make.

            evolution_precision:    The precision to require of the evolution.

            orbital_period_tolerance:    The required precision with which the
                orbital period at the present age must be reproduced.

            spin_tolerance:    The required precision with which the
                primary spin period at the present age must be reproduced.

            initial_eccentricity:    The initial eccentricity with which to
                start the evolution.

        Returns:
            None
        """

        self.binary = None
        self._best_initial_conditions = None

        if planet_formation_age:
            self.planet_formation_age = planet_formation_age
        if disk_dissipation_age is not None:
            self.disk_dissipation_age = disk_dissipation_age
        self.evolution_max_time_step = evolution_max_time_step
        self.evolution_precision = evolution_precision
        self.orbital_period_tolerance = orbital_period_tolerance
        self.spin_tolerance = spin_tolerance
        self.initial_eccentricity = initial_eccentricity

        self.target = None
        self.primary = None
        self.secondary = None

    def stellar_wsurf(self,
                      wdisk,
                      orbital_period_guess,
                      return_difference=False):
        """
        The stellar spin frquency when reproducing current porb.

        Args:
            disk_frequency:    The angular velocity of the star when it forms.

            orbital_period_guess:    A best guess value for the initial orbital
                period.

            return_difference:    If True, instead of the actual stellar
                angular velocity, the function returns the difference from the
                observed value.

        Returns:
            float or (float, float, float):

                * The angular velocity with which the star spins at the present
                  age for an evolution scenario which reproduces the current
                  orbital period. Or the difference between the spin frequency
                  and the target spin frequency if return_difference is True.

                The following are returned only if return_difference is False:

                    * The initial orbital period which reproduces the specified
                      final orbital period as close as possible.

                    * The closest final orbital period found (starting with
                      porb_initial).
        """

        disk_period = 2.0 * scipy.pi / wdisk

        porb_min, porb_max = self._find_porb_range(orbital_period_guess,
                                                   disk_period)

        if scipy.isnan(porb_min):
            assert scipy.isnan(porb_max)
            return (scipy.nan if return_difference else (scipy.nan,
                                                         scipy.nan,
                                                         scipy.nan))
        porb_initial = brentq(
            lambda porb_initial: self._try_initial_conditions(
                porb_initial,
                disk_period,
            )[0] - self.target.Porb,
            porb_min,
            porb_max,
            xtol=self.orbital_period_tolerance,
            rtol=self.orbital_period_tolerance
        )

        porb_final, spin_period = self._try_initial_conditions(
            porb_initial,
            disk_period,
        )

        spin_frequency = 2.0 * scipy.pi / spin_period

        if not return_difference:
            return spin_frequency, porb_initial, porb_final

        result = spin_frequency - 2.0 * scipy.pi / self.target.Psurf
        if (
                abs(result)
                <
                #False positive
                #pylint: disable=no-member
                abs(self._best_initial_conditions.spin_frequency
                    -
                    2.0 * scipy.pi / self.target.Psurf)
                #pylint: enable=no-member
        ):
            self._best_initial_conditions.spin_frequency = spin_frequency
            self._best_initial_conditions.orbital_period = porb_final
            self._best_initial_conditions.initial_orbital_period = porb_initial
            self._best_initial_conditions.disk_period = disk_period
        return result

    def __call__(self, target, star, planet):
        """
        Find initial conditions which reproduce the given system now.

        Args:
            target:    The target configuration to reproduce by tuning the the
                initial conditions for. The following attributes must be
                defined:

                    - age:
                        The age at which the system configuration is known.

                    - Porb:
                        The orbital period to reproduce.

                    - Psurf | Pdisk | Wdisk:
                        The stellar surface spin period to reproduce or the
                        disk locking period or the disk locking frequency.

                The following optional attributes can be specified:

                    - planet_formation_age:
                        The age at which the planet forms in Gyrs. If not
                        specified the planet is assumed to form either
                        '-past_lifetime' Gyrs before '-age' or as soon as
                        the disk dissipates.

                    - past_lifetime:
                        An alternative way of specifying when the planet
                        forms. If the '-planet_formation-age' attribute is
                        not defined, the planet is assumed to form this many
                        Gyr before '-age'.

                    - disk_dissipation_age:
                        The age at which the disk dissipates in Gyrs. If not
                        specified, it must have been defined when this solver
                        was initialized.
            star:    The star to use in the evolution, should be instance of
                evolve_interface.EvolvingStar and its dissipative properties
                should be defined.

            planet:    The planet to use in the evolution. Should be an instance
                of evolve_interface.LockedPlanet

        Returns:
            (float, float):
                * Initial orbital period.

                * Initial disk period if matching observed stellar spin or
                  stellar spin if initial disk period is specified.

            Further, the solver object has an attribute named 'binary' (an
            instance of (evolve_interface.Binary) which was evolved from
            the initial conditions found to most closely reproduce the
            specified target configuration.
        """

        def reset_best_initial_conditions():
            """Reset the entries in self._best_initial_conditions."""

            self._best_initial_conditions = Structure(
                spin_frequency=scipy.inf,
                orbital_period=scipy.nan,
                initial_orbital_period=scipy.nan,
                disk_period=scipy.nan
            )

        def get_initial_grid():
            """Tabulate stellar spin errors for a grid of disk periods."""

            reset_best_initial_conditions()

            wdisk_grid = [-200.0 * scipy.pi,
                          -20.0 * scipy.pi,
                          -2.0 * scipy.pi,
                          -0.2 * scipy.pi,
                          -0.02 * scipy.pi,
                          -0.002 * scipy.pi,
                          0.002 * scipy.pi,
                          0.02 * scipy.pi,
                          0.2 * scipy.pi,
                          2.0 * scipy.pi,
                          20.0 * scipy.pi,
                          200.0 * scipy.pi]

            stellar_wsurf_residual_grid = [
                self.stellar_wsurf(wdisk, target.Porb, True)
                for wdisk in wdisk_grid
            ]

            print('##    %25s %25s\n' % ('disk_period',
                                         'current_stellar_spin'))
            for wdisk, wsurf_residual in zip(wdisk_grid,
                                             stellar_wsurf_residual_grid):
                print('##    %25.16e %25.16e\n'
                      %
                      (2.0 * scipy.pi / wdisk,
                       wsurf_residual + 2.0 * scipy.pi / self.target.Psurf))
            print('## Target current stellar spin: '
                  +
                  repr(2.0 * scipy.pi / self.target.Psurf)
                  +
                  '\n')

            return wdisk_grid, stellar_wsurf_residual_grid

        self.target = target
        self.primary = star
        self.secondary = planet

        if not hasattr(self.target, 'disk_dissipation_age'):
            self.target.disk_dissipation_age = self.disk_dissipation_age
        if not hasattr(self.target, 'planet_formation_age'):
            if hasattr(self.target, 'past_lifetime'):
                self.target.planet_formation_age = (self.target.age
                                                    -
                                                    self.target.past_lifetime)
            else:
                self.target.planet_formation_age = self.planet_formation_age

        if not hasattr(target, 'Psurf'):
            wdisk = (target.Wdisk if hasattr(target, 'Wdisk')
                     else 2.0 * scipy.pi / target.Pdisk)
            wstar, porb_initial = self.stellar_wsurf(wdisk, target.Porb)[:2]
            return porb_initial, 2.0 * scipy.pi / wstar

        wdisk_grid, stellar_wsurf_residual_grid = get_initial_grid()

        nsolutions = 0
        for i in range(len(wdisk_grid)-1):
            if (
                    stellar_wsurf_residual_grid[i]
                    *
                    stellar_wsurf_residual_grid[i+1]
                    <
                    0
            ):
                wdisk = brentq(
                    f=self.stellar_wsurf,
                    a=wdisk_grid[i],
                    b=wdisk_grid[i+1],
                    args=(target.Porb, True),
                    xtol=self.spin_tolerance,
                    rtol=self.spin_tolerance
                )

                nsolutions += 1

            #False positive
            #pylint: disable=no-member
            return (self._best_initial_conditions.initial_orbital_period,
                    self._best_initial_conditions.disk_period)
            #pylint: enable=no-member
#pylint: enable=too-many-instance-attributes
