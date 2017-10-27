/**\file
 * 
 * \brief Declares the test suite that exercises the OrbitSolver class and
 * some other clasess necessary to accomplish this.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_ORBIT_SOLVER_H
#define __TEST_ORBIT_SOLVER_H

#include "../Evolve/OrbitSolver.h"
#include "../Star/EvolvingStar.h"
#include "../Planet/LockedPlanet.h"
#include "../Evolve/DiskBinarySystem.h"
#include "shared/PolynomialEvolution.h"
#include "shared/Common.h"
#include <iostream>
#include <vector>

namespace Evolve {

    ///Define identifiers for the quantities whose evolution we check.
    enum RealEvolutionQuantity {
        SEMIMAJOR,
        ECCENTRICITY,
        CONV_INCLINATION,
        RAD_INCLINATION,
        CONV_PERIAPSIS,
        RAD_PERIAPSIS,
        CONV_ANGMOM,
        RAD_ANGMOM,
        AGE,
        NUM_REAL_QUANTITIES
    };

    /**\brief A class that can be passed to the solution testing function
     * instead of the solver that transforms the solutions before testing.
     *
     * The tranformed orbit is obtained by applying a list of
     * OneArgumentDiffFunction objects to each of the tabulated variables.
     * Each transformation function is used only in a pre-specified age
     * range.
     *
     * \ingroup Evolve_group
     * \ingroup UnitTest_group
     */
    class TransformedSolution {
    private:
        ///The transformed orbit values.
        std::vector< const std::list<double> *> __transformed_orbit;

        ///The transformed derivatives.
//        std::vector< const std::list<double> *>  __transformed_deriv;

        ///\brief The functions to use to transform the solution.
        ///
        ///Each vector entry corresponds to one quantity in the order defined
        ///by RealEvolutionQuantity and the inner list contains the various
        ///patches to apply.
        std::vector< std::list<const Core::OneArgumentDiffFunction *> >
            __transforms;

        ///The boundaries between consecutive trasformation functions.
        std::list<double> __change_ages;

     public:
        ///\brief Default constructor, use add_transformation repeatedly to
        ///build the final solution.
        TransformedSolution() :
            __transformed_orbit(NUM_REAL_QUANTITIES, NULL),
//            __transformed_deriv(NUM_REAL_QUANTITIES)
            __transforms(AGE)
        {};

        ///Create a single piece transformed solution.
        TransformedSolution(
            ///See transforms argument to add_transformation()
            const std::vector<const Core::OneArgumentDiffFunction *> &
                transforms,

            ///The age at which this transformation starts to apply
            double start_age
        ) :
            __transformed_orbit(NUM_REAL_QUANTITIES, NULL),
//            __transformed_deriv(NUM_REAL_QUANTITIES)
            __transforms(AGE)
        {
            add_transformation(transforms, start_age);
        }

        ///Add more pieces to the transformation.
        void add_transformation(
            ///The functions to use to transform the real valued quantities.
            //The order is defined by RealEvolutionQuantity.
            const std::vector<const Core::OneArgumentDiffFunction *> &
                transforms,

            ///The age up to which this transformation applies.
            double change_age
        );

        ///\brief Apply this transformatiot to the given solution.
        ///
        ///Returns a reference to member variable, so either copy the result
        ///or do not destroy this object before use.
        const std::vector< const std::list<double> * > &operator()(
            ///Entries are assumed ordered by RealEvolutionQuantity
            const std::vector< const std::list<double> * > &solution
        );

        ///\brief The last transformed solution.
        ///
        ///Returns a reference to member variable, so either copy the result
        ///or do not destroy this object before use.
        const std::vector< const std::list<double> * > 
            &get_transformed_solution() const
            {return __transformed_orbit;}

        ///The derivative of a transformed variable at the tabulated ages.
/*        const std::list<double> *get_tabulated_var_deriv(
                EvolVarType var_type) const;*/

        ~TransformedSolution();
    };

    /**\brief Some evolution mode that changes at specified ages.
     *
     * \ingroup UnitTests_group
     */
    template<typename MODE_TYPE>
        class ExpectedEvolutionMode {
        private:
            ///The ages at which evolution mode changes occur.
            std::list<double> __age_breaks;

            ///The evolution modes that apply between consecutive __age_breaks.
            std::list<MODE_TYPE> __expected_mode;

            ///The precision with which breaks should be detected
            double __break_precision;
        public:
            ///Create.
            ExpectedEvolutionMode(double break_precision = 1e-5) :
                __break_precision(break_precision)
            {}

            ///Add an evolution mode that applies up to the given age.
            void add_break(double age, MODE_TYPE mode)
            {__age_breaks.push_back(age); __expected_mode.push_back(mode);}

            ///Is the given age close to a break (hence ambigous mode).
            bool near_break(double age) const
            {
                for(
                    std::list<double>::const_iterator 
                        break_i = __age_breaks.begin();
                    break_i != __age_breaks.end();
                    ++break_i
                )
                    if(check_diff(age, *break_i, 0.0, __break_precision))
                        return true;
                return false;
            }

            ///The evolution mode that corresponds to the given age.
            MODE_TYPE operator()(double age) const
            {
                typename std::list<MODE_TYPE>::const_iterator
                    mode_i = __expected_mode.begin();
                for(
                    std::list<double>::const_iterator age_i = __age_breaks.begin();
                        age_i != __age_breaks.end();
                        age_i++
                ) {
                    std::list<double>::const_iterator next_age_i=age_i;
                    ++next_age_i;
                    if(
                        *age_i <= age
                        &&
                        (next_age_i == __age_breaks.end() || *next_age_i > age)
                    )
                        return *mode_i;
                    ++mode_i;
                }
                std::ostringstream msg;
                msg << "Age "
                    << age
                    << " outside the range for which expected evolution modes "
                    << "are defined.";
                throw Core::Error::BadFunctionArguments(msg.str());
            }
        };

    /**\brief The test suite that exercises the OrbitSolver class.
     *
     * \ingroup UnitTests_group
     */
    class test_OrbitSolver : public Test::Suite {
    private:
        ///The solver used for the current test.
        Evolve::OrbitSolver *__solver;

        ///The system being evolved by the current test.
        Evolve::DiskBinarySystem *__system;

        ///The star used in the current test.
        Star::InterpolatedEvolutionStar *__star;

        ///Make __star a non-dissipative star with the given properties.
        void make_const_lag_star(
            const StellarEvolution::Interpolator &evolution,
            double wind_strength,
            double wind_sat_freq,
            double coupling_timescale,
            double phase_lag = 0
        );

        ///\brief Create __star with constant dissipation in a range, quickly
        ///decaying outside of that.
        void make_single_component_star(
            ///The stellar evolution to use.
            const StellarEvolution::Interpolator &evolution,

            ///The strength of the wind.
            double wind_strength,

            ///The wind saturation frequency.
            double wind_sat_freq,

            ///The core-envelope coupling timescale.
            double coupling_timescale,

            ///The minimum frequency at which the dissipation should be at its
            ///maximum value.
            double min_frequnecy,

            ///The maximum frequency at which the dissipation should be at its
            ///maximum value.
            double max_frequnecy,

            ///The scale on which frequnecy should decay.
            double decay_scale,

            ///The phase lag of the only dissipative tidal component.
            double phase_lag = 1.0e-5
        );

        StellarEvolution::MockStellarEvolution *make_no_evolution(
            double Rstar = 1.0,
            double Iconv = 1.0
        );
        StellarEvolution::MockStellarEvolution *make_linear_I_evolution();

        ///Add a planet to the given star and evolve, returning the solver.
        void evolve(
            double wdisk,
            double tdisk,
            double initial_a,
            const double *initial_Lstar,
            double initial_incl = 0.0,
            double planet_mass = 1.0,
            ///If NaN defaults to tdisk.
            double tplanet = Core::NaN,
            double max_age = MAX_AGE,
            double planet_radius = 1.0,
            double precision = 1e-6,
            double max_step_factor = 1e-3
        );

        ///Return the last calculated evolution.
        std::vector< const std::list<double> *> get_evolution() const;

        ///\brief Tests the latest evolution calculated by the solver against
        ///the given tracks.
        void test_solution(
            const std::vector< const std::list<double> * > &
                tabulated_real_quantities,
            std::vector<const Core::OneArgumentDiffFunction *>
                expected_real_quantities,
            const ExpectedEvolutionMode<Core::EvolModeType> &
                expected_evol_mode,
            const ExpectedEvolutionMode<bool> &expected_wind_mode,
            double min_age,
            double max_age,
            bool debug_mode = false
        );

        ///\brief Test a planet-less scenario computed in 3 different ways:
        ///1) withou a planet 2) without dissipation 3) with massless planet
        void test_no_planet_scenario(
            const StellarEvolution::Interpolator &stellar_evol,
            double *initial_Lstar,
            double windK,
            double wind_sat_freq,
            double core_env_coupling_time,
            std::vector<const Core::OneArgumentDiffFunction *>
                &expected_evolution,
            const ExpectedEvolutionMode<bool> &expected_wind_mode,
            double max_age = MAX_AGE,
            bool debug_mode = false
        );
    protected:
        ///No fixtures at this time
        void setup() {};

        ///No fixtures at this time
        void tear_down() {};
    public:
        ///Create the test suite.
        test_OrbitSolver();

        ///\brief Tests the evolution of stellar systems with disk locking
        ///times longer than the stellar lifetimes for a non-evolving star.
        void test_disk_locked_no_stellar_evolution();

        ///\brief Tests the evolution of stellar systems with disk locking
        ///times longer than the stellar lifetimes, for evolving stars.
        void test_disk_locked_with_stellar_evolution();

        ///\brief Tests the rotational evolution of a star without a planet
        ///and without a disk locking phase.
        ///
        ///Calculates each evolution with NO_PLANET evolution mode, and with
        ///FAST_PLANET/SLOW_PLANET evolution mode, but zero planet mass.
        void test_no_planet_evolution();

        ///\brief Tests the evolution of the orbit plus stellar rotation,
        ///starting with the planet already present and ensuring it does not
        ///die or lock to the star.
        void test_unlocked_evolution();

        ///\brief Tests the evolution of the orbit plus stellar rotation for
        ///the case where the star is locked in synchronous rotation with the
        ///orbit.
        void test_locked_evolution();

        ///\brief Tests an evolution that starts locked to a disk then
        ///immediately is locked to the planet, which then falls below the
        ///Roche radius.
        void test_disklocked_to_locked_to_noplanet();

        ///\brief Tests an evolution that starts locked to a disk then
        ///the orbital period is shorter than the stellar spin period and
        ///eventually the planet falls below the Roche radius.
        void test_disklocked_to_fast_to_noplanet();

        ///\brief Tests an evolution that starts locked to a disk then
        ///the orbital period is shorter than the stellar spin period and
        ///then is locked.
        void test_disklocked_to_fast_to_locked();

        ///\brief Tests an evolution that starts locked to a disk then
        ///immediately is locked to the planet and then the lock is broken
        ///leaving the orbital period shorter than the stellar spin period.
        void test_disklocked_to_locked_to_fast();

        ///\brief Tests an evolution with only the 1-0 component in a polar
        ///orbit.
        ///
        ///Expectations are that no evolution will occur.
        void test_polar_1_0_evolution();

        ///\brief Tests an evolution with only the 2-0 component starting
        ///in a polar orbit.
        ///
        ///Approximately the expected evolution has the star spinning down
        ///linearly with time and the inclination should follow from angular
        ///momentum conservation with the orbit approximately not changing.
        void test_polar_2_0_evolution();

        ///\brief Test an evolution with only the 1-0 component but in an
        ///arbitrarily inclined orbit.
        ///
        ///From angular momentum conservation, the obliquity should satisfy:
        /// \f$ \cos\theta  = \frac{T^2 - S^2 - L^2}{2 S L} \f$
        ///Where:
        /// - S: stellar spin angular momentum
        /// - L: orbital angular momentum
        /// - T: total angular momentum (stellar and orbital)
        ///
        ///From Lai 2012:
        ///  - \f$ \dot{S} = -\sigma c^2 (1 - c^2) \f$
        ///  - \f$ \dot{c} = \frac{\sigma}{S} c^2 (1 - c^2)\left[c + \frac{S}{L} \right] \f$
        ///
        ///Per mathematica, the evolution of S is given by:
        /// \f$ \frac{1}{2} L^3 \left(-\frac{4}{L^2+S^2-T^2}-\frac{2 \log \left(L^2+S^2-T^2\right)}{L^2-T^2}-\frac{-\frac{2 L^3 \log \left(-L^2+S^2+T^2\right)}{L^2-T^2}+(L+T) \log (-L+S-T)+(L-T) (\log (L+S-T)+\log (-L+S+T))+(L+T) \log (L+S+T)}{L T^2}\right) \f$
        void test_oblique_1_0_evolution();
    };//End test_OrbitSolver class.

}//End Evolve namespace.


#endif
