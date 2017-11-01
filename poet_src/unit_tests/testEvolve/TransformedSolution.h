/**\file
 *
 * \brief Declare a class that transforms an orbital solution before checking.
 *
 * This is useful if an analytical solution exists for some function of the
 * evolved quantities, which cannot be analytically inverted.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TRANSFORMED_SOLUTION_H
#define __TRANSFORMED_SOLUTION_H

#include "RealEvolutionQuantity.h"
#include "../../Core/Functions.h"
#include <vector>
#include <list>

namespace Evolve {

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

} //End Evolve namespace.

#endif
