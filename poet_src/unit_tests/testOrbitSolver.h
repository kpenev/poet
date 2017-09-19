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

#if 0
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
        std::vector< std::list<double> > __transformed_orbit,

            ///The transformed derivatives.
            __transformed_deriv;

        ///The functions to use to transform the semimajor axis.
        std::list<const OneArgumentDiffFunction *> __a_transforms,

            ///\brief The function to use  to transfom the convective zone
            ///angular momentum.
            __Lconv_transforms,
            
            ///\brief The function to use  to transfom the radiative zone
            ///angular momentum.
            __Lrad_transforms;

        ///The ages at which the orbit was tabulated.
        const std::list<double> *__ages;

        ///The boundaries between consecutive trasformation functions.
        std::list<double> __change_ages;

        ///The evolution mode at each tabulated age
        const std::list<EvolModeType> *__evolution_mode;

        ///Applies a list of transformations to a variable.
        void transform(
                ///The variable to transform.
                const std::list<double>* var,

                ///The derivative of the variable to transform.
                const std::list<double>* deriv,
                
                ///The list of trasformation functions.
                const std::list<const OneArgumentDiffFunction*> &transforms,
                
                ///Which variable is being transformed.
                EvolVarType var_type);
    public:
        ///\brief Default constructor, use add_transformation repeatedly to
        ///build the final solution.
        TransformedSolution() :
            __transformed_orbit(3),
            __transformed_deriv(3)
        {};

        ///Create a single piece transformed solution.
        TransformedSolution(
            ///The function to use to transform the semimajor axis.
            const OneArgumentDiffFunction &a_transform,

            ///The function to use to transform the convective zone
            ///angular momentum.
            const OneArgumentDiffFunction &Lconv_transform,

            ///The function to use to transform the radiative zone
            ///angular momentum.
            const OneArgumentDiffFunction &Lrad_transform,

            ///The age at which this transformation starts to apply
            double start_age
        ) :
            __transformed_orbit(3),
            __transformed_deriv(3)
        {
            add_transformation(a_transform,
                               Lconv_transform,
                               Lrad_transform,
                               start_age);
        }

        ///Add more pieces to the transformation.
        void add_transformation(
            ///The function to use to transform the semimajor axis after the
            ///last change age up to the given one.
            const OneArgumentDiffFunction &a_transform,

            ///The function to use to transform the convective zone angular
            ///momentum after the last change age up to the given one.
            const OneArgumentDiffFunction &Lconv_transform,

            ///The function to use to transform the radiative zone angular
            ///momentum after the last change age up to the given one.
            const OneArgumentDiffFunction &Lrad_transform,

            ///The age up to which this transformation applies.
            double change_age
        );

        ///\brief The solver to apply this transformation to.
        ///
        ///It should already have a solution stored.
        void operator()(const OrbitSolver &solver);

        ///The value of a transformed variable at the tabulated ages.
        const std::list<double>
            *get_tabulated_var(
                    ///Which variable to return.
                    EvolVarType var_type) const;

        ///The derivative of a transformed variable at the tabulated ages.
        const std::list<double> *get_tabulated_var_deriv(
                EvolVarType var_type) const;

        ///Returns a list of the evolution modes.
        const std::list<EvolModeType> *get_tabulated_evolution_mode() const
        {return __evolution_mode;} 
    };
#endif

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
                        (next_age_i == __age_breaks.end() || *next_age_i >= age)
                    )
                        return *mode_i;
                    ++mode_i;
                }
                assert(false);
            }
        };

    /**\brief The test suite that exercises the OrbitSolver class.
     *
     * \ingroup UnitTests_group
     */
    class test_OrbitSolver : public Test::Suite {
    private:
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
            NUM_REAL_EVOL_QUANTITIES
        };

        ///The solver used for the current test.
        Evolve::OrbitSolver *__solver;

        ///The system being evolved by the current test.
        Evolve::DiskBinarySystem *__system;

        ///The star used in the current test.
        Star::InterpolatedEvolutionStar *__star;

        friend std::ostream &operator<<(std::ostream &os,
                                        RealEvolutionQuantity q);

        ///Make __star a non-dissipative star with the given properties.
        void make_const_lag_star(
            const StellarEvolution::Interpolator &evolution,
            double wind_strength,
            double wind_sat_freq,
            double coupling_timescale,
            double phase_lag = 0
        );

        StellarEvolution::MockStellarEvolution *make_no_evolution();
        StellarEvolution::MockStellarEvolution *make_linear_I_evolution();

        ///Add a planet to the given star and evolve, returning the solver.
        void evolve(
            double wdisk,
            double tdisk,
            double initial_a,
            double *initial_Lstar,
            double initial_incl = 0.0,
            double planet_mass = 1.0,
            ///If NaN defaults to tdisk.
            double tplanet = Core::NaN,
            double max_age = MAX_AGE
        );

        ///\brief Tests the latest evolution calculated by the solver against
        ///the given tracks.
        void test_solution(
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
//        void test_unlocked_evolution();

        ///\brief Tests the evolution of the orbit plus stellar rotation for
        ///the case where the star is locked in synchronous rotation with the
        ///orbit.
//        void test_locked_evolution();

        ///\brief Tests an evolution that starts locked to a disk then
        ///immediately is locked to the planet, which then falls below the
        ///Roche radius.
//        void test_disklocked_to_locked_to_noplanet();

        ///\brief Tests an evolution that starts locked to a disk then
        ///the orbital period is shorter than the stellar spin period and
        ///eventually the planet falls below the Roche radius.
//        void test_disklocked_to_fast_to_noplanet();

        ///\brief Tests an evolution that starts locked to a disk then
        ///the orbital period is shorter than the stellar spin period and
        ///then is locked.
//        void test_disklocked_to_fast_to_locked();

        ///\brief Tests an evolution that starts locked to a disk then
        ///immediately is locked to the planet and then the lock is broken
        ///leaving the orbital period shorter than the stellar spin period.
//        void test_disklocked_to_locked_to_fast();
    };//End test_OrbitSolver class.

}//End Evolve namespace.


#endif
