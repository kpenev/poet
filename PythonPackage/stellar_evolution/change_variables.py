"""Allows finding mass and age given other properties at fixed [Fe/H]."""
import sys
sys.path.append('..')

from stellar_evolution.library_interface import MESAInterpolator
from stellar_evolution.derived_stellar_quantities import\
    TeffK,\
    LogGCGS,\
    RhoCGS
from basic_utils import Structure
import scipy
import scipy.interpolate
import scipy.linalg
import scipy.optimize

class QuantityEvaluator :
    """Evaluate stellar quantities fixing out of range issues etc."""

    def __init__(self, interpolator, metallicity, min_mass = 0.0, max_mass
                 = scipy.inf, **reference_values) :
        """Set-up to use the given interpolator."""

        self.interp = interpolator
        self.metallicity = metallicity
        self.reference_values = reference_values
        self.min_mass = min_mass
        self.max_mass = max_mass
        for var in ['teff', 'logg', 'lum', 'rho'] :
            if var not in reference_values :
                self.reference_values[var] = 0.0

    def _evaluate(self, quantity, age) :
        """Evaluate the given quantity at the given age."""

        if age < quantity.min_age or age > quantity.max_age :
            return scipy.nan
        else :
            return quantity(age)

    def teff(self, mass, age) :
        """Return the effective temperature for the given stellar params."""

        if mass < self.min_mass or mass > self.max_mass : return scipy.nan
        return self._evaluate(
            TeffK(self.interp('radius', mass, self.metallicity), 
                  self.interp('lum', mass, self.metallicity)),
            age
        ) - self.reference_values['teff']

    def rho(self, mass, age) :
        """Return the density for the given stellar params."""

        if mass < self.min_mass or mass > self.max_mass : return scipy.nan
        return self._evaluate(
            RhoCGS(mass, self.interp('radius', mass, self.metallicity)), 
            age
        ) - self.reference_values['rho']

    def logg(self, mass, age) :
        """Return log10(surface gravity) for the given stellar params."""

        if mass < self.min_mass or mass > self.max_mass : return scipy.nan
        return self._evaluate(
            LogGCGS(mass, self.interp('radius', mass, self.metallicity)), 
            age
        ) - self.reference_values['logg']

    def lum(self, mass, age) :
        """Return log10(surface gravity) for the given stellar params."""

        if mass < self.min_mass or mass > self.max_mass : return scipy.nan
        return self._evaluate(
            self.interp('lum', mass, self.metallicity), 
            age
        ) - self.reference_values['lum']

class VarChangingInterpolator(MESAInterpolator) :
    """
    Enhance interpolators to find mass and age given other properties.

    Instances have the following attributes:
        - grid:
            A structure with attributes: 
                - masses:
                    Stellar masses of the grid nodes at which the dependent
                    variables are known.
                - ages:
                    Stellar ages of the grid nodes at which the dependent
                    variables are known.
                - metallicities:
                    Metallicities of the grid nodes at which the dependent
                    variables are known.
                - weights:
                    The weights to use when interpolating the grid to a
                    specified metallicity.
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
            mass, followed by age and finally metallicity.

            The grid must be fine enough to ensure that no grid cell entirely
            contains an iso-contour of either dependent variable where
            interpolation will be attempted.
    """

    def _get_quantity(self, name, mass, metallicity) :
        """
        Return a quantity at the given mass and metallicity.

        Args:
            - name:
                The name of the quantity to return.
            - mass:
                The stellar mass for which this quantity should apply
            - metallicity:
                The [Fe/H] for which this quantity should apply.

        Returns:
            A callable returning the value of the quantity at a given age.
            Also has min_age and max_age attributes defining the range over
            which it is defined.
        """

        if name == 'teff' :
            return TeffK(self('radius', mass, metallicity),
                         self('lum', mass, metallicity))
        elif name == 'logg' :
            return LogGCGS(mass,
                           self('radius', mass, metallicity))
        elif name == 'rho' :
            return RhoCGS(mass,
                          self('radius', mass, metallicity))
        else :
            return self(name, mass, metallicity)

    def search_near(self, mass, age, metallicity, **kwargs) :
        """
        Search for mass & age near the given ones to match two other vars.

        Args:
            - mass:
                The value of the mass to search near.
            - age:
                The value of the age to search near.
            - metallicity:
                The value of [Fe/H] at which the variable change is taking
                place.

        Keyword only arguments, must be exactly two of the following:
            - teff:
                The effective temperature to match.
            - logg:
                The log10(gravitation acceleration) to match.
            - lum:
                The luminosity to match.
            - rho:
                The density to match.

        Returns:
            - (mass, age):
                The mass and age where the given dependent variables are
                matched.
        """

        assert(len(kwargs) == 2)
        evaluator = QuantityEvaluator(self,
                                      metallicity,
                                      min_mass = self.grid.masses[0],
                                      max_mass = self.grid.masses[-1],
                                      **kwargs)
        missmatch = [getattr(evaluator, quantity_name)
                     for quantity_name, reference in kwargs.items()]
        solution = scipy.optimize.root(
            lambda m_t: scipy.array([miss(*m_t) for miss in missmatch]),
            [mass, age],
            method = 'lm',
            options = dict(ftol = 1e-15, xtol = 1e-15)
        )
        if solution.success : return solution.x
        else : return None

    def _add_grid_variable(self, variable) :
        """
        Adds another dependent variable to self.grid.

        Args:
            - variable:
                The name of the variable to add. See class documentation for
                details.

        Returns: None
        """

        setattr(
            self.grid,
            variable,
            scipy.zeros((self.grid.masses.size,
                         self.grid.ages.size,
                         self.grid.metallicities.size))
        )

        for metallicity_index, metallicity in enumerate(
                self.grid.metallicities
        ) :
            for mass_index, mass in enumerate(self.grid.masses) :
                quantity = self._get_quantity(variable, mass, metallicity)
                age_in_range = scipy.logical_and(
                    self.grid.ages > quantity.min_age,
                    self.grid.ages < quantity.max_age
                )
                getattr(self.grid, variable)[
                    mass_index,
                    age_in_range,
                    metallicity_index
                ] = quantity(self.grid.ages[age_in_range])
                if not self._defined_weights :
                    self.grid.weights[mass_index,
                                      :,
                                      metallicity_index] = age_in_range
                else : 
                    assert(
                        (self.grid.weights[mass_index, :, metallicity_index]
                         ==
                         age_in_range).all()
                    )

    def _interpolate_grid_variable(self, var_name, metallicity) :
        """
        Interpolate one of the grid variables to a specified metallicity.

        Args:
            - var_name:
                The name of the variable to interpolate.
            - metallicity:
                The [Fe/H] value to inteprolate to.

        Returns: A 2-D scipy array contaning the interepolated variable at
        the grid masses and ages. 
        """

        result = scipy.empty((self.grid.masses.size, self.grid.ages.size))
        for mass_index in range(self.grid.masses.size) :
            for age_index in range(self.grid.ages.size) :
                interp_y = getattr(self.grid, var_name)[mass_index,
                                                        age_index,
                                                        :]
                weights = self.grid.weights[mass_index, age_index, :]
                if weights.sum() < 2 : 
                    result[mass_index, age_index] = scipy.nan
                else :
                    result[
                        mass_index,
                        age_index
                    ] = scipy.interpolate.InterpolatedUnivariateSpline(
                        x = self.grid.metallicities,
                        y = interp_y,
                        k = 1
                    )(metallicity)
        return result

    def __init__(self, **kwargs) :
        """
        Prepare an interpolator able to find mass, age from other quantities.

        Args:
            - change_from:
                A length-2 iterable listing the two quantities which will
                be transformed to mass and age. Any two of the following:
                'teff', 'logg', 'lum',  'rho'. This is case-insensitive.

        Keyword only arguments: see MESAInterpolator.__init__

        Returns: None
        """

        super().__init__(**kwargs)
        self.grid = Structure(
            metallicities = scipy.linspace(
                float(self.track_metallicities[0]),
                float(self.track_metallicities[-1]),
                3 * len(self.track_metallicities)
            ),
            masses = scipy.linspace(
                float(self.track_masses[0]),
                float(self.track_masses[-1]),
                10 * len(self.track_masses)
            ),
            ages = scipy.linspace(1e-2, 13.71, 412),
        )
        self.grid.weights = scipy.empty((self.grid.masses.size,
                                         self.grid.ages.size,
                                         self.grid.metallicities.size),
                                        dtype = bool)
        self._defined_weights = False

    def change_variables(self, metallicity, **kwargs) :
        """
        Change from two of (lum, rho, logg, Teff) to mass & age.

        Args:
            - metallicity:
                The value of [Fe/H] at which this variable change is taking
                place.

        Keyword only arguments, must be exactly two of the following:
            - teff:
                The effective temperature to match.
            - logg:
                The log10(gravitation acceleration) to match.
            - lum:
                The luminosity to match.
            - rho:
                The density to match.

        Returns: 
            The mass and age at which the keyword arguments are matched
            as list of tuplse of (mass, age).
        """

        def find_candidate_cells() :
            """
            Identify grid cells possibly containing a solution.

            Args: None

            Returns:
                - possible_solutions:
                    A 2-D array matching the shape of grid variables with
                    True entries for grid cells which may contain a solution.
                - var_diff:
                    A 3-D array with the last two dimensions matching the
                    shape of grid variables, consisting of two slabs giving
                    the difference from the target values of each of the
                    grid variables used for variable change.
            """

            possible_solutions = True
            var_diff = scipy.empty(
                (2, self.grid.masses.size, self.grid.ages.size)
            )
            for var_index, (var_name, var_value) in enumerate(
                    kwargs.items()
            ) :
                if not hasattr(self.grid, var_name) :
                    self._add_grid_variable(var_name)
                var_diff[var_index] = self._interpolate_grid_variable(
                    var_name,
                    metallicity
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
                sign_change = scipy.logical_or(
                    scipy.logical_or(sign_change_right[:, :-1],
                                     sign_change_up[:-1, :]),
                    scipy.logical_or(sign_change_right[:, 1:],
                                     sign_change_up[1:, :])
                )
                possible_solutions = scipy.logical_and(possible_solutions,
                                                       sign_change)

            return possible_solutions, var_diff

        def bilinear_coef_equations(m0, m1, t0, t1) :
            """
            Equations for the coef. of a bi-linear func. over mass/age cell.

            Args:
                - m0:
                    The lower mass boundary of the cell.
                - m1:
                    The upper mass boundary of the cell.
                - t0:
                    The lower age boundary of the cell.
                - t1:
                    The upper age boundary of the cell.

            Returns: A 2-D scipy array containing the matrix defining the
            equations for the coefficients of a bi-linear function over the
            specified cell.
            """

            coef_equations = scipy.empty((4, 4))
            coef_equations[:, 0] = 1
            coef_equations[0:2, 1] = m0
            coef_equations[2:4, 1] = m1
            coef_equations[0::2, 2] = t0
            coef_equations[1::2, 2] = t1
            coef_equations[0, 3] = m0 * t0
            coef_equations[1, 3] = m0 * t1
            coef_equations[2, 3] = m1 * t0
            coef_equations[3, 3] = m1 * t1
            return coef_equations

        def find_roots(coef) :
            """
            Return the simultaneous roots of two bilinear functions.

            Args:
                - coef: 
                    The coefficients of the two bilinear functions. Should be
                    a 2-D scipy array with the outer index iterating over the
                    function and the inner indices iterating over the
                    coefficients of the corresponding bilinear function.
            
            Returns: A list of 2-tuples contaning the simultaneous zeros
            of the two functions.
            """

            a = scipy.linalg.det(coef[:,2:])
            b = (scipy.linalg.det(coef[:,::3])
                 -
                 scipy.linalg.det(coef[:,1:3])) 
            c = scipy.linalg.det(coef[:,:2])
            det = b * b - 4.0 * a * c
            if det < 0 : return []
            sqrt_det = det**0.5
            time_roots = (scipy.array([(-b - sqrt_det), (-b + sqrt_det)])
                          /
                          (2.0 * a))
            mass_roots =  - (
                (coef[0, 0] + coef[0, 2] * time_roots)
                /
                (coef[0, 1] + coef[0, 3] * time_roots)
            )
            return zip(mass_roots, time_roots)

        def get_unique_solutions(solutions, telorance = 1e-8) :
            """Merge repeated solutions from the given list."""

            ref_ind = 0
            while ref_ind < len(solutions) - 1 :
                ref_solution = solutions[ref_ind]
                target_ind = ref_ind + 1
                while target_ind < len(solutions) :
                    target_solution = solutions[target_ind]
                    if (
                            scipy.linalg.norm(
                                (ref_solution - target_solution)
                                /
                                (ref_solution + target_solution)
                            ) 
                            < 
                            tolerance
                    ) :
                        del solutions[target_ind]

        possible_solutions, var_diff = find_candidate_cells()
        result = []
        for mass_index, age_index in scipy.transpose(
                scipy.nonzero(possible_solutions)
        ) :
            t0, t1 = self.grid.ages[age_index: age_index + 2]
            m0, m1 = self.grid.masses[mass_index: mass_index + 2]
            if scipy.isnan(
                    var_diff[:,
                             mass_index : mass_index + 2,
                             age_index : age_index + 2]
            ).any() :
                continue
            coef_equations = bilinear_coef_equations(m0, m1, t0, t1)
            coef = scipy.empty((2, 4))
            for var_index in range(2) :
                coef[var_index] = scipy.linalg.solve(
                    coef_equations, 
                    var_diff[var_index,
                             mass_index : mass_index + 2,
                             age_index : age_index + 2].flatten()
                )
            candidate_solutions = find_roots(coef)
            for mass, age in candidate_solutions :
                if (
                        mass >= m0 and mass <= m1 
                        and
                        age >= t0 and age <= t1
                ) :
                    solution = self.search_near(mass = mass,
                                                age = age,
                                                metallicity = metallicity,
                                                **kwargs)
                    if solution is not None : 
                        result.append(solution)
        return result
