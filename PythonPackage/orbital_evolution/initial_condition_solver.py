from math import pi

class InitialConditionSolver :
    """Find initial conditions which reproduce a given system now."""

    def _try_initial_conditions(self, initial_orbital_period, disk_period) :
        """
        Get present orbital and stellar spin periods for initial conditions.

        Args:
            - initial_orbital_period:
                The initial orbital period to calculate the deviation for.
            - disk_period:
                The disk locking period to calculate the deviation for.

        Returns:
            - orbital_period:
                The present day orbital period of the system resulting when
                the evolution is started with the input periods.
            - spin_period:
                The present day surface spin of the star resulting when the
                evolution is started with the input periods.
        """

        self.binary.primary.select_interpolation_region(
            self.disk_dissipation_age
        )
        self.binary.primary.detect_stellar_wind_saturation()

        self.binary.configure(self.primary.core_formation_age,
                              float('nan'),
                              float('nan'),
                              numpy.array([0.0]),
                              None,
                              None,
                              'LOCKED_SURFACE_SPIN')
        self.binary.secondary.configure(
            self.planet_formation_age,
            self.binary.primary.mass,
            self.semimajor(initial_orbital_period),
            0.0,
            numpy.array([0.0]),
            None,
            None,
            False,
            True,
            True
        )
        self.binary.evolve(
            self.system_parametrs.age,
            self.evolution_max_time_step,
            self.evolution_precision,
            None
        )
        final_state = self.binary.final_state()
        assert(final_state.age == self.system_parametrs.age)
        return (self.binary.orbital_period(final_state.semimajor),
                (2.0 * pi
                 *
                 self.binary.primary.envelope_inertia(final_state.age)
                 /
                 final_state.envelope_angmom))

    def _find_porb_range(self, guess_porb_initial, disk_period) :
        """
        Find initial orbital period range where final porb error flips sign.

        Args:
            - guess_porb_initial:
                An initial guess for where the sign change occurs.

            - disk_period:
                The disk locking period to assume during the search.

        Returns:
            A pair of initial orbital periods for which the sign of the final
            orbital period error changes.
        """

        porb_min, porb_max = numpy.nan, numpy.nan
        porb_initial = guess_porb_initial
        porb, psurf = self._try_initial_conditions(porb_initial, disk_period)
        porb_error = porb - self.system.Porb
        guess_porb_error = porb_error
        step = 2.0 if guess_porb_error < 0 else 0.5

        while porb_error * guess_porb_error > 0 and porb_initial < 100.0 :
            if porb_error < 0 : porb_min = porb_initial
            else : porb_max = porb_initial
            porb_initial *= step
            porb, psurf = self.__try_initial_conditions(porb_initial,
                                                        disk_period)
            if not numpy.isnan(porb) :
                porb_error = porb - self.system.Porb

        if numpy.isnan(porb_error) : return numpy.nan, numpy.nan

        if porb_error < 0 : porb_min = porb_initial
        else :
            porb_max = porb_initial
            if porb_error == 0 : porb_min = porb_initial

        return porb_min, porb_max

    def __init__(self,
                 past_lifetime = None,
                 timeout_seconds = None) :
        """
        Initialize the object.

        Args:
            - past_lifetime: If not None, the system is started at an
                             age = current age - past_lifetime, and the
                             initial conditions at that age are solved for.
                             Otherwise, the starting age is left to the poet
                             command line arguments.
            - timeout_seconds: If an orbital evolution takes longer than this
                               many seconds to finish it is canceled and
                               invalid results are assumed. Use None to
                               disable timing out.

        Returns: None.
        """

        self.__past_lifetime = past_lifetime
        self.timeout_seconds = timeout_seconds

    def stellar_wsurf(self,
                      wdisk,
                      orbital_period_guess,
                      return_difference=False) :
        """
        The stellar spin frquency when reproducing current porb.

        Args:
            - disk_frequency: The angular velocity of the star when it forms.
            - orbital_period_guess: A best guess value for the initial
                                    orbital period.
            - return_difference: If True, instead of the actual stellar
                                 angulare velocity, the function returns the
                                 difference from the observed value.

        Returns: The angular velocity with which the star spins at the
                 present age for an evolution scenario which reproduces the
                 current orbital period.
        """

        disk_period = 2.0 * pi / wdisk

        porb_min, porb_max = self._find_porb_range(orbital_period_guess,
                                                   disk_period)

        if numpy.isnan(porb_min) :
            assert(numpy.isnan(porb_max))
            return (numpy.nan if return_difference else (numpy.nan,
                                                         numpy.nan,
                                                         numpy.nan))
        porb_initial = brentq(
            lambda porb_initial : self.__try_initial_conditions(
                porb_initial,
                disk_period,
            )[0] - self.system.Porb,
            porb_min,
            porb_max
        )

        orbital_period, spin_period = self.__try_initial_conditions(
            porb_initial,
            disk_period,
        )

        spin_frequency = 2.0 * pi / spin_period

        if not return_difference : 
            return spin_frequency, porb_initial, orbital_period

        result = spin_frequency - 2.0 * pi / self.system.Psurf
        if (
                abs(result) 
                < 
                abs(self.__best_initial_conditions.spin_frequency 
                    - 
                    2.0 * pi / self.system.Psurf) 
        ) :
            self.__best_initial_conditions.spin_frequency = spin_frequency
            self.__best_initial_conditions.orbital_period = orbital_period
            self.__best_initial_conditions.initial_orbital_period = (
                porb_initial
            )
            self.__best_initial_conditions.disk_period = disk_period
        return result

    def __call__(self, system) :
        """
        Find initial conditions which reproduce the given system now.

        Args:
            - system: The system to find the initial conditions for. Should
                      be generated by calling one of :
                        - hatsouth_to_solver_system()
                        - exoplanets_org_solver_system().
                        - nasa_solver_system()

        Returns:
            - porb: Initial orbital period.
            - pdisk: Initial disk period.
        """

        def reset_best_initial_conditions() :
            """Reset the entries in self.__best_initial_conditions."""

            self.__best_initial_conditions=Structure(
                spin_frequency = numpy.inf,
                orbital_period = numpy.nan,
                initial_orbital_period = numpy.nan,
                disk_period = numpy.nan
            )

        def get_initial_grid() :
            """Tabulate stellar spin errors for a grid of disk periods."""

            reset_best_initial_conditions()

            wdisk_grid = [-200.0*pi,
                          -20.0*pi,
                          -2.0*pi,
                          -0.2*pi,
                          -0.02*pi,
                          -0.002*pi,
                          0.002*pi,
                          0.02*pi,
                          0.2*pi,
                          2.0*pi,
                          20.0*pi,
                          200.0*pi]

            stellar_wsurf_residual_grid = [
                self.stellar_wsurf(wdisk, system.Porb, True)
                for wdisk in wdisk_grid
            ]

            self.__scenario_file.write(
                '##    %25s %25s\n'%('disk_period', 'current_stellar_spin')
            )
            for wdisk, wsurf_residual in zip(wdisk_grid,
                                             stellar_wsurf_residual_grid) :
                self.__scenario_file.write(
                    '##    %25.16e %25.16e\n'
                    %
                    (2.0 * pi / wdisk, 
                     wsurf_residual + 2.0 * pi / self.system.Psurf)
                )
            self.__scenario_file.write('## Target current stellar spin: '
                                       +
                                       repr(2.0 * pi / self.system.Psurf)
                                       +
                                       '\n')

            return wdisk_grid, stellar_wsurf_residual_grid

        self.system = system
        self.planet_formation_age = (
            getattr(self.system_parameters, 'start_age', None) 
            if self._past_lifetime is None else
            (self.system_parametrs.age - self._past_lifetime)
        )

        wdisk_grid, stellar_wsurf_residual_grid = get_initial_grid()

        nsolutions = 0
        for i in range(len(wdisk_grid)-1) :
            if (
                    stellar_wsurf_residual_grid[i]
                    *
                    stellar_wsurf_residual_grid[i+1]
                    <
                    0
            ) :
                wdisk = brentq(
                    f = self.stellar_wsurf,
                    a = wdisk_grid[i],
                    b = wdisk_grid[i+1],
                    args = (system.Porb,
                            True),
                    xtol = 1e-8,
                    rtol = 1e-8
                )

                nsolutions += 1
