#!/usr/bin/python3 -u

from matplotlib import pyplot
import sys
sys.path.append('../PythonPackage')

from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import *
from stellar_evolution.change_variables import QuantityEvaluator
from math import pi
import numpy

def get_solutions(interp,
                  target_mass,
                  target_metallicity,
                  target_age) :
    """
    Return all solutions found attempting to recover the specified params.

    Args:
        - interp:
            The interpolator on which variable changes are based. 
        - target_mass:
            The true stellar mass of the correct solution.
        - target_metallicity:
            The true stellar metallicity of the correct solution.
        - target_age:
            The true age of the correct solution.

    Returns: A possibly empty list of 2-tuples of (mass, age) recovered.
    """

def visualize_mismatch(interp,
                       target_mass,
                       target_metallicity,
                       target_age) :
    """
    Visualize a solution which does not match expectations.

    For now assumes that the input variables are stellar density and
    effective temperature.

    Args:
        - interp:
            The interpolator on which variable changes are based. 
        - target_mass:
            The true stellar mass of the correct solution.
        - target_metallicity:
            The true stellar metallicity of the correct solution.
        - target_age:
            The true age of the correct solution.

    Returns: None
    """

    target_radius = interp('radius', target_mass, target_metallicity)
    target_lum = interp('lum', target_mass, target_metallicity)
    target_teff = TeffK(target_radius, target_lum)
    target_rho = RhoCGS(target_mass, target_radius)
    solutions = interp.change_variables(
        metallicity = target_metallicity,
        rho = target_rho(target_age),
        teff = target_teff(target_age)
    )

    for solution_mass, solution_age in solutions :
        radius = interp('radius', solution_mass, target_metallicity)
        lum = interp('lum', solution_mass, target_metallicity)
        teff = TeffK(radius, lum)
        rho = RhoCGS(solution_mass, radius)
        print(
            'M* = %s, t = %s, Rho* = %s (%s), Teff = %s (%s)'
            %
            (
                repr(solution_mass),
                repr(solution_age),
                repr(rho(solution_age)),
                repr(target_rho(target_age)),
                repr(teff(solution_age)),
                repr(target_teff(target_age))
            )
        )
        mass_diff = target_mass - solution_mass
        age_diff = target_age - solution_age
        plot_x = numpy.linspace(-2.0, 2.0, 1000)
        mass_age_track = [
            (target_mass + x * mass_diff, target_age + x * age_diff)
            for x in plot_x
        ]
        evaluator = QuantityEvaluator(interp, target_metallicity)
        plot_teff = numpy.array([
            evaluator.teff(m, t) - target_teff(target_age)
            for m, t in mass_age_track
        ])
        plot_rho = numpy.array([
            evaluator.rho(m, t) - target_rho(target_age)
            for m, t in mass_age_track
        ])

        figure, axis1 = pyplot.subplots()
        axis1.plot(plot_x, plot_teff, label = 'Teff')
        axis1_yrange = max([abs(y) for y in axis1.get_ylim()])
        axis1.set_ylim(-axis1_yrange, axis1_yrange)
        axis2 = axis1.twinx()
        axis2.plot(plot_x, plot_rho, label = 'rho')
        axis2_yrange = max([abs(y) for y in axis2.get_ylim()])
        axis2.set_ylim(-axis2_yrange, axis2_yrange)
        pyplot.show()
        pyplot.cla()

        pyplot.plot(plot_x, plot_teff**2 + plot_rho**2)
        pyplot.show()
        pyplot.cla()

        m = numpy.linspace(1.5 * solution_mass - 0.5 * target_mass,
                           2.5 * target_mass - 1.5 * solution_mass,
                           100)
        t = numpy.linspace(1.5 * solution_age - 0.5 * target_age,
                           2.5 * target_age - 1.5 * solution_age,
                           100)
        grid_m, grid_t = numpy.meshgrid(m, t)
        square_diff = numpy.empty((m.size, t.size))
        teff_ref = target_teff(target_age)
        rho_ref = target_rho(target_age)
        for m_ind, m_val in enumerate(m) :
            for t_ind, t_val in enumerate(t) :
                square_diff[m_ind, t_ind] = (
                    (evaluator.teff(m_val, t_val) - teff_ref)**2
                    +
                    (evaluator.rho(m_val, t_val) - rho_ref)**2
                )
        pyplot.imshow(square_diff,
                      interpolation = 'bilinear',
                      origin = 'lower',
                      aspect = 'equal')
        pyplot.show()

def random_tests(interp) :
    """Run a series of randomly selected tests."""

    line_fmt = '%25s: %25s %25s %25s %25s %25s %25s'
    print(line_fmt % ('[Fe/H]',
                      'M*(target)',
                      't(target)',
                      'DM*',
                      'dt',
                      'M*(found)',
                      't(found)'))
    for target_mass in numpy.random.uniform(0.4, 1.2, 10) :
        for target_metallicity in numpy.random.uniform(-1.0, 0.5, 10) :
            target_radius = interp('radius', target_mass, target_metallicity)
            target_lum = interp('lum', target_mass, target_metallicity)
            target_teff = TeffK(target_radius, target_lum)
            target_rho = RhoCGS(target_mass, target_radius)
            for target_age in numpy.random.uniform(target_radius.min_age,
                                                   min(target_radius.max_age,
                                                       13.7),
                                                   10) :
                solutions = interp.change_variables(
                    metallicity = target_metallicity,
                    rho = target_rho(target_age),
                    teff = target_teff(target_age)
                )
                if not solutions : 
                    found_mass_str = found_age_str = mass_diff_str = \
                    age_diff_str = '-'
                else :
                    min_distance = numpy.inf
                    for solution_mass, solution_age in solutions :
                        distance = ((solution_mass - target_mass)**2
                                    +
                                    (solution_age - target_age)**2)
                        if distance < min_distance :
                            found_mass_str = repr(solution_mass)
                            found_age_str = repr(solution_age)
                            mass_diff_str = repr(target_mass - solution_mass)
                            age_diff_str = repr(target_age - solution_age)
                            min_distance = distance
                print(line_fmt
                      %
                      (repr(target_metallicity),
                       repr(target_mass),
                       repr(target_age),
                       mass_diff_str,
                       age_diff_str,
                       found_mass_str,
                       found_age_str))

if __name__ == '__main__' :
    manager = StellarEvolutionManager('../stellar_evolution_interpolators')
    interp = manager.get_interpolator_by_name('default')
    visualize_mismatch(interp = interp,
                       target_mass = 0.77285588027585062, 
                       target_metallicity = -0.28278013953588377,
                       target_age = 13.812435550152532)
#    random_tests(interp)
