#!/usr/bin/python3 -u

from matplotlib import pyplot
import sys
sys.path.append('../PythonPackage')

from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.derived_stellar_quantities import *
from stellar_evolution.change_variables import QuantityEvaluator
from math import pi
import numpy

if __name__ == '__main__' :
    manager = StellarEvolutionManager('../stellar_evolution_interpolators')
    interp = manager.get_interpolator_by_name('default')

    target_mass = 0.53
    target_metallicity = 0.0
    target_age = pi
    target_radius = interp('radius', target_mass, target_metallicity)
    target_lum = interp('lum', target_mass, target_metallicity)
    target_teff = TeffK(target_radius, target_lum)
    target_rho = RhoCGS(target_mass, target_radius)

    solutions = interp.change_variables(metallicity = target_metallicity,
                                        rho = target_rho(target_age),
                                        teff = target_teff(target_age))
    for mass, age in solutions :
        radius = interp('radius', mass, target_metallicity)
        lum = interp('lum', mass, target_metallicity)
        teff = TeffK(radius, lum)
        rho = RhoCGS(mass, radius)
        print(
            'M* = %s, t = %s, Rho* = %s (%s), Teff = %s (%s)'
            %
            (
                repr(mass),
                repr(age),
                repr(rho(age)),
                repr(target_rho(target_age)),
                repr(teff(age)),
                repr(target_teff(target_age))
            )
        )
        mass_diff, age_diff = target_mass - mass, target_age - age
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

        m = numpy.linspace(1.5 * mass - 0.5 * target_mass,
                           2.5 * target_mass - 1.5 * mass,
                           100)
        t = numpy.linspace(1.5 * age - 0.5 * target_age,
                           2.5 * target_age - 1.5 * age,
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
