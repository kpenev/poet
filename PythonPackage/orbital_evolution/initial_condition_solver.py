
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
         MAKE IT POSSIBLE TO GET FINAL CONFIG <++>

        return orbital_period, spin_period

    def __find_porb_range(self, guess_porb_initial, disk_period) :
        """Find a range of initial porb where final porb error flips sign."""

        porb_min, porb_max = numpy.nan, numpy.nan
        porb_initial=guess_porb_initial
        porb, psurf = self.__try_initial_conditions(porb_initial,
                                                    disk_period)
        porb_error = porb - self.system.Porb
        guess_porb_error = porb_error
        step = 2.0 if guess_porb_error < 0 else 0.5

        while porb_error*guess_porb_error > 0 and porb_initial < 100.0 :
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
                 poet_cmdline,
                 scenario_file,
                 input_columns,
                 output_columns,
                 trial_evolutions_fname,
                 selected_evolution_fname,
                 past_lifetime = None,
                 timeout_seconds = None) :
        """
        Initialize the object.

        Args:
            - poet_cmdline: A command line with which to start the POET
                            process.
            - scenario_file: An open file object which poet is reading for
                             scenarios to try.
            - input_columns: A list of the columns expected by poeit in the
                             scenario file.
            - trial_evolutions_fname: A %-substitution pattern for the
                                      filenames to use for each trial
                                      evolution while searching for the
                                      solution. The pattern may contain a
                                      %(trial) substition, which gets filled
                                      by an integer incremented at each trial
                                      step. If not, the file gets overwritten
                                      at each trial. These files are NOT
                                      automatically deleted.
            - selected_evolution_fname: The filename to save the best
                                        evolution found as. It should not
                                        coincide with any possible resolution
                                        of trial_evolutions_fname.
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

        self.__poet_cmdline = poet_cmdline
        self.__scenario_file = open(scenario_file, 'w')
        self.__scenario_file.write(
            "# PID %d POET command line: '" % current_process().pid
            + 
            "' '".join(poet_cmdline)
            +
            "'\n"
        )
        self.__scenario_file.flush()

        self.poet = Popen(poet_cmdline, stdin=PIPE, bufsize=0)
        self.__input_columns = input_columns
        self.__output_columns = output_columns
        self.__trial_index = 0
        self.__past_lifetime = past_lifetime
        self.trial_evolutions_fname = trial_evolutions_fname
        self.selected_evolution_fname = selected_evolution_fname
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

        porb_min, porb_max = self.__find_porb_range(orbital_period_guess,
                                                    disk_period)

        self.__scenario_file.write(
            '# For disk period = %25.16e initial orbital period range is '
            '[%25.16e, %25.16e].\n'
            %
            (disk_period, porb_min, porb_max)
        )
        self.__scenario_file.flush()

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

        self.__scenario_file.write(
            '# For disk period = %25.16e the current orbital period of '
            '%25.16e is reproduced by starting from orbital period of '
            '%25.16e. Current stellar angular velocity = %25.16e.\n'
            %
            (disk_period, self.system.Porb, porb_initial, spin_frequency)
        )

        if not return_difference : 
            shutil.copyfile(self.__evol_fname, self.selected_evolution_fname)
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
            shutil.copyfile(self.__evol_fname, self.selected_evolution_fname)
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

        try :
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

                    self.__scenario_file.write(
                        '###    Solution %d: \n'
                        '###        age = %25.16e\n'
                        '###        P0 = %25.16e\n'
                        '###        pdisk = %25.16e\n'
                        '###        Porb = %25.16e\n'
                        '###        Wsurf = %25.16e\n'
                        %
                        (nsolutions,
                         self.system.age,
                         self.__best_initial_conditions.\
                            initial_orbital_period,
                         self.__best_initial_conditions.disk_period,
                         self.__best_initial_conditions.orbital_period,
                         self.__best_initial_conditions.spin_frequency)
                    )
                    nsolutions += 1
            if nsolutions == 0 :
                    self.__scenario_file.write(
                        '###    No solution found: \n'
                        '###        age = %25.16e\n'
                        '###        P0 = %25.16e\n'
                        '###        pdisk = %25.16e\n'
                        '###        Porb = %25.16e\n'
                        '###        Wsurf = %25.16e\n'
                        %
                        (self.system.age,
                         self.__best_initial_conditions.\
                            initial_orbital_period,
                         self.__best_initial_conditions.disk_period,
                         self.__best_initial_conditions.orbital_period,
                         self.__best_initial_conditions.spin_frequency)
                    )
        except :
            self.__scenario_file.write(
                '#' + format_exc().replace('\n', '\n#')
            )
            raise
