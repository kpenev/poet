"""Define class using machine learning for faster initial condition solving."""

import numpy

from orbital_evolution import \
    InitialConditionSolver,\
    EvolvingStar,\
    LockedPlanet,\
    lgQ

class LearningICSolver(InitialConditionSolver):
    """Use machine learning wrapped around I.C. solver to speed up solution."""

    #pylint: disable=too-many-branches
    def _get_porb_solver_ml_features(self, disk_period, final_porb=None):
        """
        Return numpy array of values required by the ML to guess Porb range.

        Args:
            disk_period(float):    The disk period the evolution will be
                starting with.

            final_porb(float):    The final orbital period that must be
                reproduced. By default ``self.target.Porb`` is used.

        Returns:
            numpy.array:
                The values of the features machine learning uses to predict the
                initial orbital period range. The first value is always the
                final orbital period.
        """

        feh = None
        result = [final_porb or self.target.Porb, self.target.age]
        for body in [self.primary, self.secondary]:
            for _, dissipation in sorted(body.dissipation.items()):
                result.append(lgQ(dissipation['reference_phase_lag']))
                if dissipation.tidal_frequency_breaks is not None:
                    result.extend(dissipation['tidal_frequency_breaks'])
                result.extend(dissipation['tidal_frequency_powers'])
                if dissipation.spin_frequency_breaks is not None:
                    result.extend(dissipation['spin_frequency_breaks'])
                result.extend(dissipation['spin_frequency_powers'])
                if isinstance(body, LockedPlanet):
                    result.append(body.radius)

            result.append(body.mass)
            if isinstance(body, EvolvingStar):
                if 'metallicity' not in self._exclude_ml_param:
                    if feh is None:
                        result.append(body.metallicity)
                    else:
                        assert body.metallicity == feh
                for param in ['wind_strength',
                              'wind saturation frequency',
                              'diff rot coupling timescale']:
                    if param not in self._exclude_ml_param:
                        result.append(getattr(body, param.replace(' ', '_')))

        for param, value in sorted(self.parameters.items()):
            if param not in self._exclude_ml_param:
                if param == 'initial secondary angmom':
                    result.extend(value)
                else:
                    result.append(value)

        result.append(disk_period)
        return numpy.array(result)
    #pylint: enable=too-many-branches


    def _try_initial_conditions(self,
                                initial_orbital_period,
                                disk_period,
                                save=False):
        """
        Pass training data to machine learning model for initial conditions.

        See `InitialConditionSolver._try_initial_conditions()` for full
        documentation.
        """

        final_porb, final_primary_pspin = super()._try_initial_conditions(
            initial_orbital_period,
            disk_period,
            save
        )
        if final_porb:
            ml_features = self._get_porb_solver_ml_features(disk_period,
                                                            final_porb)
            if self.validate_training_data(ml_features, initial_orbital_period):
                self._guess_porb_range.add_training_data(
                    ml_features,
                    initial_orbital_period
                )
        return final_porb, final_primary_pspin


    def _find_porb_range(self, guess_porb_initial, disk_period):
        """
        ML aided finding of period range where final porb error flips sign.

        See `InitialConditionSolver._find_porb_range()` for full description.

        If the machine learning produces estimated period range,
        ``guess_porb_initial`` is ignored.
        """

        ml_porb_range = self._guess_porb_range(
            self._get_porb_solver_ml_features(disk_period)
        )
        if ml_porb_range[0] is None:
            return super()._find_porb_range(guess_porb_initial, disk_period)

        assert ml_porb_range[1] > ml_porb_range[0]
        period_search_factor = self.period_search_factor
        self.period_search_factor = ml_porb_range[1] / ml_porb_range[0]
        result = super()._find_porb_range(ml_porb_range[0], disk_period)
        self.period_search_factor = period_search_factor
        return result


    def __init__(self, *, exclude_ml_param, **kwargs):
        """
        Prepare the solver.

        Args:
            exclude_ml_param(iterable):    List of properties not to include in
                the features that are used for machine learning. Should be a
                subset of:
                    'metallicity',
                    'wind strength',
                    'wind saturation frequency',
                    'diff rot coupling timescale',
                    'initial inclination',
                    'initial eccentricity',
                    'initial secondary angmom',
                    'planet formation age',
                    'disk dissipation age'
                Any values not in the above list are ignored.

            kwargs:    Passed directly to `InitialConditionSolver.__init__()`

        Returns:
            None
        """

        self._exclude_ml_param = exclude_ml_param
        super().__init__(**kwargs)
        #TODO: replace with ML
        self._guess_porb_range = lambda _: (None, None)


    #Meant to be overwritten if needed
    #pylint: disable=no-self-use
    #pylint: disable=unused-argument
    def validate_training_data(self, ml_features, initial_porb):
        """Return True iff the given training data should be used."""

        return True
    #pylint: enable=no-self-use
    #pylint: enable=unused-argument
