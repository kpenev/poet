"""Allows finding mass and age given other properties at fixed [Fe/H]."""

import sys

import numpy
import scipy.interpolate
import scipy.linalg
import scipy.optimize

sys.path.append('..')

#Need to add POET package path to module search path before importing
#pylint: disable=wrong-import-position
from stellar_evolution.library_interface import MESAInterpolator
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS
from basic_utils import Structure
#pylint: enable=wrong-import-position

class QuantityEvaluator:
    """Evaluate stellar quantities fixing out of range issues etc."""

    def __init__(self,
                 interpolator,
                 feh=0.0,
                 **reference_values):
        """Set-up to use the given interpolator."""

        self.interp = interpolator
        self.feh = feh
        self.reference_values = reference_values
        for var in ['teff', 'logg', 'lum', 'rho']:
            if var not in reference_values:
                self.reference_values[var] = 0.0

    @staticmethod
    def _evaluate(quantity, age):
        """Evaluate the given quantity at the given age."""

        if quantity.min_age < age < quantity.max_age:
            return quantity(age)
        return numpy.nan

    def teff(self, mass, age, feh=None):
        """Return the effective temperature for the given stellar params."""

        if feh is None:
            feh = self.feh

        if self.interp.in_range(mass, feh):
            radius = self.interp('radius', mass, feh)
            luminosity = self.interp('lum', mass, feh)
            result = self._evaluate(
                TeffK(radius, luminosity),
                age
            ) - self.reference_values['teff']
            radius.delete()
            luminosity.delete()
            return result

        return numpy.nan

    def rho(self, mass, age, feh=None):
        """Return the density for the given stellar params."""

        if feh is None:
            feh = self.feh

        if self.interp.in_range(mass, feh):
            radius = self.interp('radius', mass, self.feh)
            result = self._evaluate(
                RhoCGS(mass, radius),
                age
            ) - self.reference_values['rho']
            radius.delete()
            return result

        return numpy.nan

    def logg(self, mass, age, feh=None):
        """Return log10(surface gravity) for the given stellar params."""

        if feh is None:
            feh = self.feh

        if self.interp.in_range(mass, feh):
            radius = self.interp('radius', mass, self.feh)
            result = self._evaluate(
                LogGCGS(mass, radius),
                age
            ) - self.reference_values['logg']
            radius.delete()
            return result

        return numpy.nan

    def lum(self, mass, age, feh=None):
        """Return log10(surface gravity) for the given stellar params."""

        if feh is None:
            feh = self.feh

        if self.interp.in_range(mass, feh):
            luminosity = self.interp('lum', mass, self.feh)
            result = self._evaluate(
                luminosity,
                age
            ) - self.reference_values['lum']
            luminosity.delete()
            return result

        return numpy.nan

class VarChangingInterpolator(MESAInterpolator):
    """
    Enhance interpolators to find mass and age given other properties.

    Attributes:
        grid:    A structure with attributes:
            - masses:
                Stellar masses of the grid nodes at which the dependent
                variables are known.

            - ages:
                Stellar ages of the grid nodes at which the dependent
                variables are known.

            - feh:
                [Fe/H] values of the grid nodes at which the dependent
                variables are known.

            - weights:
                The weights to use when interpolating the grid to a
                specified [Fe/H].

            - teff (may not be present):
                The stellar effective temperature in Kelvin at the grid
                nodes.

            - logg (may not be present):
                The log10(stellar gravity in cgs) at the grid nodes.

            - lum (may not be present):
                The stellar luminosity, in solar luminosities at the grid
                nodes.

            - rho(may not be present):
                The stellar density in cgs at the grid nodes.

            The first index for each variable or the weights attribute is
            mass, followed by age and finally [Fe/H].

    Notes:
        The grid must be fine enough to ensure that no grid cell entirely
        contains an iso-contour of either dependent variable where
        interpolation will be attempted.
    """

    def _get_quantity(self, name, mass, feh):
        """
        Return a quantity at the given mass and [Fe/H].

        Args:
            name:    The name of the quantity to return.

            mass:    The stellar mass for which this quantity should apply

            feh:    The [Fe/H] for which this quantity should apply.

        Returns:
            callable:
                A callable returning the value of the quantity at a given age.
                Also has min_age and max_age attributes defining the range over
                which it is defined.
        """

        if name == 'teff':
            return TeffK(self('radius', mass, feh),
                         self('lum', mass, feh))

        if name == 'logg':
            return LogGCGS(mass,
                           self('radius', mass, feh))

        if name == 'rho':
            return RhoCGS(mass,
                          self('radius', mass, feh))

        return self(name, mass, feh)

    def search_near(self, mass, age, feh, **kwargs):
        """
        Search for mass & age near the given ones to match two other vars.

        Args:
            mass:    The value of the mass to search near.

            age:    The value of the age to search near.

            feh:    The value of [Fe/H] at which the variable change is taking
                place.

            kwargs:    must be exactly two of the following:

                - teff:
                    The effective temperature to match.

                - logg:
                    The log10(gravitation acceleration) to match.

                - lum:
                    The luminosity to match.

                - rho:
                    The density to match.

        Returns:
            (float, float):
                The mass and age where the given dependent variables are
                matched.
        """

        assert len(kwargs) == 2
        #False positive, members created by setattr, so pylint does not see them
        #pylint: disable=no-member
        evaluator = QuantityEvaluator(self,
                                      feh,
                                      min_mass=self.grid.masses[0],
                                      max_mass=self.grid.masses[-1],
                                      **kwargs)
        #pylint: enable=no-member
        missmatch = [getattr(evaluator, quantity_name)
                     for quantity_name, reference in kwargs.items()]
        solution = scipy.optimize.root(
            lambda m_t: numpy.array([miss(*m_t) for miss in missmatch]),
            [mass, age],
            method='lm',
            options=dict(ftol=1e-15, xtol=1e-15)
        )
        if solution.success:
            return solution.x

        return None

    def _add_grid_variable(self, variable):
        """
        Adds another dependent variable to self.grid.

        Args:
            variable:    The name of the variable to add. See class
                documentation for details.

        Returns:
            None
        """

        #False positive, members created by setattr, so pylint does not see them
        #pylint: disable=no-member
        setattr(
            self.grid,
            variable,
            numpy.zeros((self.grid.masses.size,
                         self.grid.ages.size,
                         self.grid.feh.size))
        )

        for feh_index, feh in enumerate(self.grid.feh):
            for mass_index, mass in enumerate(self.grid.masses):
                quantity = self._get_quantity(variable, mass, feh)
                age_in_range = numpy.logical_and(
                    self.grid.ages > quantity.min_age,
                    self.grid.ages < quantity.max_age
                )
                getattr(self.grid, variable)[
                    mass_index,
                    age_in_range,
                    feh_index
                ] = quantity(self.grid.ages[age_in_range])

                if not self.defined_weights:
                    self.grid.weights[mass_index,
                                      :,
                                      feh_index] = age_in_range
                else:
                    assert(
                        (self.grid.weights[mass_index, :, feh_index]
                         ==
                         age_in_range).all()
                    )

        self.defined_weights = True
        #pylint: enable=no-member

    def _interpolate_grid_variable(self, var_name, feh):
        """
        Interpolate one of the grid variables to a specified [Fe/H].

        Args:
            var_name:    The name of the variable to interpolate.

            feh:    The [Fe/H] value to inteprolate to.

        Returns:
            2-D numpy array:
                The interepolated variable at the grid masses and ages.
        """

        #False positive, members created by setattr, so pylint does not see them
        #pylint: disable=no-member
        result = numpy.empty((self.grid.masses.size, self.grid.ages.size))
        for mass_index in range(self.grid.masses.size):
            for age_index in range(self.grid.ages.size):
                interp_y = getattr(self.grid, var_name)[mass_index,
                                                        age_index,
                                                        :]
                weights = self.grid.weights[mass_index, age_index, :]
                if weights.sum() < 2:
                    result[mass_index, age_index] = numpy.nan
                else:
                    result[
                        mass_index,
                        age_index
                    ] = scipy.interpolate.InterpolatedUnivariateSpline(
                        x=self.grid.feh,
                        y=interp_y,
                        k=1
                    )(feh)
        return result
        #pylint: enable=no-member

    def _define_var_change_grid(self, feh, masses, ages):
        """
        Create a new grid with the given locations of the nodes.

        Creates self.grid with all arguments as same-name members and an
        additional weight member containing an empty array to fill with weights
        later when grid variables start to be calculated. Also creates
        self._defined_weigths to keep track if weights have been previously
        initialized.

        Args:
            feh:    The [Fe/H] values at which to tabulate the dependent
                variables.

            masses:    The stellar masses at which to tabulate the dependent
                variables.

            - ages:    The ages (in Gyrs) at which to tabulate the dependent
                variables.

        Returns:
            None
        """

        self.grid = Structure(feh=feh,
                              masses=masses,
                              ages=ages)
        #False positive, members created by setattr, so pylint does not see them
        #pylint: disable=no-member
        self.grid.weights = numpy.empty(
            (
                self.grid.masses.size,
                self.grid.ages.size,
                self.grid.feh.size
            ),
            dtype=bool
        )
        #pylint: enable=no-member
        self.defined_weights = False

    def _find_candidate_cells(self, feh, **kwargs):
        """
        Identify grid cells possibly containing a solution.

        Args:
            feh:    See change_variables().

            kwargs:    See change_variables().

        Returns:
            2-D array:
                Matching the shape of grid variables with True entries for
                grid cells which may contain a solution.

            3-D array:
                The last two dimensions matching the shape of grid
                variables, consisting of two slabs giving the difference
                from the target values of each of the grid variables used
                for variable change.
        """

        possible_solutions = True

        #False positive, members created by setattr, so pylint does not
        #see them
        #pylint: disable=no-member
        var_diff = numpy.empty(
            (2, self.grid.masses.size, self.grid.ages.size)
        )
        #pylint: enable=no-member
        for var_index, (var_name, var_value) in enumerate(
                kwargs.items()
        ):
            if not hasattr(self.grid, var_name):
                self._add_grid_variable(var_name)
            var_diff[var_index] = self._interpolate_grid_variable(
                var_name,
                feh
            ) - var_value
            sign_change_right = (
                (var_diff[var_index][1:, :]
                 *
                 var_diff[var_index][:-1, :])
                <=
                0
            )
            sign_change_up = (
                (var_diff[var_index][:, 1:]
                 *
                 var_diff[var_index][:, :-1])
                <=
                0
            )
            sign_change = numpy.logical_or(
                numpy.logical_or(sign_change_right[:, :-1],
                                 sign_change_up[:-1, :]),
                numpy.logical_or(sign_change_right[:, 1:],
                                 sign_change_up[1:, :])
            )
            possible_solutions = numpy.logical_and(possible_solutions,
                                                   sign_change)

        return possible_solutions, var_diff

    @staticmethod
    def _bilinear_coef_equations(m_low, m_high, age_low, age_high):
        """
        Equations for the coef. of a bi-linear func. over mass/age cell.

        Args:
            m_low:    The lower mass boundary of the cell.

            m_high:    The upper mass boundary of the cell.

            age_low:    The lower age boundary of the cell.

            age_high:    The upper age boundary of the cell.

        Returns:
            2-D numpy array:
                The matrix defining the equations for the coefficients of a
                bi-linear function over the specified cell.
        """

        coef_equations = numpy.empty((4, 4))
        coef_equations[:, 0] = 1
        coef_equations[0:2, 1] = m_low
        coef_equations[2:4, 1] = m_high
        coef_equations[0::2, 2] = age_low
        coef_equations[1::2, 2] = age_high
        coef_equations[0, 3] = m_low * age_low
        coef_equations[1, 3] = m_low * age_high
        coef_equations[2, 3] = m_high * age_low
        coef_equations[3, 3] = m_high * age_high
        return coef_equations

    @staticmethod
    def _find_bilinear_roots(coef):
        """
        Return the simultaneous roots of two bilinear functions.

        Args:
            coef:    The coefficients of the two bilinear functions. Should
                be a 2-D numpy array with the outer index iterating over the
                function and the inner indices iterating over the
                coefficients of the corresponding bilinear function.

        Returns:
            [(float, float), ...]:
                A list of 2-tuples contaning the simultaneous zeros of the
                two functions.
        """

        #These are standard names for quadratic equation.
        #pylint: disable=invalid-name
        a = scipy.linalg.det(coef[:, 2:])
        b = (scipy.linalg.det(coef[:, ::3])
             -
             scipy.linalg.det(coef[:, 1:3]))
        c = scipy.linalg.det(coef[:, :2])
        #pylint: enable=invalid-name
        det = b * b - 4.0 * a * c

        if det < 0:
            return []

        sqrt_det = det**0.5
        time_roots = (numpy.array([(-b - sqrt_det), (-b + sqrt_det)])
                      /
                      (2.0 * a))
        mass_roots = -(
            (coef[0, 0] + coef[0, 2] * time_roots)
            /
            (coef[0, 1] + coef[0, 3] * time_roots)
        )
        return zip(mass_roots, time_roots)

    def __init__(self,
                 grid_feh,
                 grid_masses,
                 grid_ages,
                 **kwargs):
        """
        Prepare an interpolator able to find mass, age from other quantities.

        Keyword only arguments: see MESAInterpolator.__init__

        Returns: None
        """

        super().__init__(**kwargs)
        self.defined_weights = False
        self._define_var_change_grid(
            feh=grid_feh,
            masses=grid_masses,
            ages=grid_ages
        )

    #Attempting to simplify furthe results in a mess.
    #pylint: disable=too-many-locals
    def change_variables(self, feh, **kwargs):
        """
        Change from two of (lum, rho, logg, Teff) to mass & age.

        Args:
            feh:    The value of [Fe/H] at which this variable change is taking
                place.

            kwargs:    must be exactly two of the following:

                - teff:
                    The effective temperature to match.

                - logg:
                    The log10(gravitation acceleration) to match.

                - lum:
                    The luminosity to match.

                - rho:
                    The density to match.

        Returns:
            [(float, float), ...]:
                The mass and age at which the keyword arguments are matched
                as list of tuples of (mass, age).
        """

        possible_solutions, var_diff = self._find_candidate_cells(feh, **kwargs)
        result = []
        for mass_index, age_index in numpy.transpose(
                numpy.nonzero(possible_solutions)
        ):
            #False positive, members created by setattr, so pylint does not see
            #them
            #pylint: disable=no-member
            age_low, age_high = self.grid.ages[age_index: age_index + 2]
            m_low, m_high = self.grid.masses[mass_index: mass_index + 2]
            #pylint: enable=no-member
            if numpy.isnan(
                    var_diff[
                        :,
                        mass_index : mass_index + 2,
                        age_index : age_index + 2
                    ]
            ).any():
                continue

            coef_equations = self._bilinear_coef_equations(m_low,
                                                           m_high,
                                                           age_low,
                                                           age_high)
            coef = numpy.empty((2, 4))
            for var_index in range(2):
                coef[var_index] = scipy.linalg.solve(
                    coef_equations,
                    var_diff[var_index,
                             mass_index : mass_index + 2,
                             age_index : age_index + 2].flatten()
                )
            candidate_solutions = self._find_bilinear_roots(coef)
            for mass, age in candidate_solutions:
                if m_low <= mass <= m_high and age_low <= age <= age_high:
                    solution = self.search_near(mass=mass,
                                                age=age,
                                                feh=feh,
                                                **kwargs)
                    if solution is not None:
                        result.append(solution)
        return result
    #pylint: enable=too-many-locals
