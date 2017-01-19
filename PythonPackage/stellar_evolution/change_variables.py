"""Allows finding mass and age given other properties at fixed [Fe/H]."""

from .library_interface import MESAInterpolator
import scipy
import scipy.interpolate
import scipy.linalg

class VarChangingInterpolator(MESAInterpolator) :
    """Enhance interpolators to find mass and age given other properties."""

    def _change(grid, metallicity, **kwargs) :
        """
        Change from two dependent variables (L, rho, logg, Teff) to mass/age.

        Args:
            - grid:
                A structure with attributes named metallicities, masses and
                ages containing the metallicity, mass and age locations of
                the grid nodes at which the dependent variables are known
                (see the keyword arguments) and two dependent variables with
                the same names as the keyword arguments used contaning the
                corresponding values at the grid points. The grid must be
                fine enough to ensure that no grid cell entirely contains an
                iso-contour of either dependent variable.
            - metallicity:
                The value of [Fe/H] at which this variable change is taking
                place.

        Keyword only arguments, must be exactly two of the following:
            - teff:
                The effective temperature to match.
            - logg:
                The log10(gravitation acceleration) to match.
            - L:
                The luminosity to match.
            - rho:
                The density to match.

        Returns: 
            The mass and age at which the keyword arguments are matched. If
            multiple solutions are found the returned value is a list of
            tuples.
        """

        possible_solutions = True
        var_diff = numpy.empty((2, grid.masses.size, grid.ages.size))
        for var_index, (var_name, var_value) in enumerate(kwargs.items()) :
            var_diff[var_index] = scipy.interpolate.interp1d(
                grid.metallicities,
                getattr(grid, var_name),
                axis = 0,
                kind = 'cubic'
            )(metallicity) - var_value
            sign_change_right = var_diff[1:, :] * var_diff[:-1, :] < 0
            sign_change_up = var_diff[:, 1:] * var_diff[:, :-1] < 0
            sign_change = scipy.logical_or(
                scipy.logical_or(sign_change_right[:, :-1],
                                 sign_change_up[:-1, :]),
                scipy.logical_or(sign_change_right[:, 1:],
                                 sign_change_up[1:, :])
            )
            possible_solutions = scipy.logical_and(possible_solutions,
                                                   sign_change)
        for mass_index, age_index in scipy.transpose(
                scipy.nonzero(possible_solutions)
        ) :
            t0, t1 = grid.ages[age_index: age_index + 2]
            m0, m1 = grid.masses[mass_index: mass_index + 2]
            coef_equations = numpy.empty((4, 4))
            coef_equations[:, 0] = 1
            coef_equations[0:2, 1] = m0
            coef_equations[2:4, 1] = m1
            coef_equations[0::2, 2] = t0
            coef_equations[1::2, 2] = t1
            coef_equations[0, 3] = m0 * t0
            coef_equations[1, 3] = m1 * t0
            coef_equations[2, 3] = m0 * t1
            coef_equations[3, 3] = m1 * t1
            coef = numpy.empty((2, 4))
            for var_index in range(2) :
                coef[var_index] = scipy.linalg.solve(
                    coef_equations, 
                    var_diff[var_index,
                             mass_index : mass_index + 2,
                             age_index : age_index + 2].flatten()
                )
            a = scipy.linalg.det(coef[:,2:])
            b = scipy.linalg.det(coef[:,::3]) - scipy.linalg.det(coef[:,1:3]) 
            c = scipy.linalg.det(coef[:,:2])
