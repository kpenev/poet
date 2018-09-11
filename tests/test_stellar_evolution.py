#!/usr/bin/python3 -u

"""Unit tests for stellar evolution interpolation."""

import os.path
import unittest

import scipy
import scipy.integrate

from add_poet_to_path import poet_root
from stellar_evolution.manager import StellarEvolutionManager
from stellar_evolution.library_interface import MESAInterpolator

class TestStellarEvolution(unittest.TestCase):
    """Exercise the stellar evolution interpolation."""

    #Simple enough
    #pylint: disable=too-many-locals
    def _test_quantity_derivative(self,
                                  quantity_name,
                                  mass,
                                  feh,
                                  min_reference_age):
        """Numerically test the derivative of a quantity."""

        quantity = self.interp(quantity_name,
                               mass,
                               feh)

        check_ages = scipy.random.uniform(min_reference_age,
                                          quantity.max_age,
                                          self.num_test_ages)
        check_against = quantity.deriv(check_ages)
        reference_age = None
        for age_index, age in enumerate(check_ages):
            age_range = quantity.continuous_range(age)
            print('Age range: ' + repr(age_range))
            reference_age = max(min_reference_age,
                                age_range[0])

            initial_deriv = quantity.deriv(reference_age)
            print(quantity_name + '(%f) = %s' % (reference_age,
                                                 initial_deriv))

            for check_order in [0, 1]:
                print('Checking order %d derivative @ t = %f'
                      %
                      (check_order, age))
                expected = (
                    check_against[check_order][age_index]
                ) - initial_deriv[check_order]

                #Pylint false positive: perfectly legitimate use!
                #pylint: disable=cell-var-from-loop
                got, error = scipy.integrate.quad(
                    lambda lnt: (
                        quantity.deriv(
                            scipy.exp(lnt)
                        )[check_order + 1]
                        *
                        scipy.exp(lnt)
                    ),
                    scipy.log(reference_age),
                    scipy.log(age),
                    limit=1000
                )
                #pylint: enable=cell-var-from-loop
                error = max(3.0 * error, 1e-8 * abs(expected))
                print('Expected %f, got %f +- %f'
                      %
                      (expected, got, error))
                with self.subTest(quantity=quantity_name,
                                  M=mass,
                                  FeH=feh,
                                  t=age,
                                  order=check_order):
                    self.assertLessEqual(
                        abs(expected - got),
                        error
                    )
        quantity.delete()
    #pylint: enable=too-many-locals

    def setUp(self):
        """Read a serialized interpolator in self.interp."""

        self.interp = StellarEvolutionManager(
            os.path.join(poet_root, 'stellar_evolution_interpolators')
        ).get_interpolator_by_name('default')
        self.num_test_masses = 10
        self.num_test_feh = 10
        self.num_test_ages = 10
        print('Done')

    def test_derivatives(self):
        """Tests that derivatives integrate to lower order ones."""

        check_masses = scipy.random.uniform(
            min(self.interp.track_masses),
            max(self.interp.track_masses),
            self.num_test_masses
        )
        check_masses[0] = 1.0
        for mass in check_masses:
            check_feh = scipy.random.uniform(
                min(self.interp.track_feh),
                max(self.interp.track_feh),
                self.num_test_feh
            )
            check_feh[0] = 0.0
            for feh in check_feh:
                min_reference_age = self.interp('RADIUS',
                                                mass,
                                                feh).min_age
                for quantity_name in MESAInterpolator.quantity_list:
                    print('Checking %s(M = %f, [Fe/H] = %f)'
                          %
                          (quantity_name, mass, feh))
                    self._test_quantity_derivative(quantity_name,
                                                   mass,
                                                   feh,
                                                   min_reference_age)

if __name__ == '__main__':
    unittest.main()
