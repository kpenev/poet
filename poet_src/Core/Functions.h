/**\file
 *
 * \brief A hierarchy of classes representing functions.
 *
 * \ingroup Utilities_group
 */

#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H

#include "../Core/SharedLibraryExportMacros.h"
#include "Common.h"
#include "Error.h"
//#include <solve_polynomial.h>
#include <cmath>
#include <limits>
#include <valarray>
#include <list>
#include <iterator>
#include <iostream>

#ifndef NO_SERIALIZE
	#include <boost/serialization/base_object.hpp>
	#include <boost/archive/text_oarchive.hpp>
	#include <boost/archive/text_iarchive.hpp>
	#include <boost/serialization/export.hpp>
#endif

namespace Core {

    class InterpSolutionIterator;

    ///\brief The base class for functions which take a single argument and
    ///return a single value.
    template<class InputDataType, class OutputDataType>
        class LIB_LOCAL OneArgumentFunction {
#ifndef NO_SERIALIZE
            ///Needed for serialization to work.
            friend class boost::serialization::access;

            ///Serialize this function.
            template<class Archive>
                void serialize(Archive &, const unsigned int) {}
#endif
        public:
            ///The value of the function at the given abscissa.
            virtual OutputDataType operator()(
                InputDataType in_value
            ) const =0;

            ///The lower end of the range over which the function is defined
            virtual InputDataType range_high() const=0;

            ///The upper end of the range over which the function is defined
            virtual InputDataType range_low() const=0;

            ///\brief An iterator over the abscissas where the function takes
            ///the given y value.
            virtual InterpSolutionIterator crossings(double y=0) const =0;

            ///Provide a virtual destructor for a virtual class.
            virtual ~OneArgumentFunction() {};
        };

    ///A class representing arbitrary order derivatives of a function.
    class LIB_LOCAL FunctionDerivatives {
    public:
        ///\brief Derivative of the given order of the function with respect
        ///to its argument.
        virtual double order(unsigned deriv_order=1) const =0;

        ///Clean up.
        virtual ~FunctionDerivatives() {};
    };

    ///A class for the derivatives of a cubic spline (=0 for order>2).
    class LIB_LOCAL CubicSplineDerivatives : public FunctionDerivatives {
    private:
        double zeroth, ///< The value of the function
               first, ///< The first derivative
               second; ///< The second derivative
    public:
        ///\brief Constuct a spline derivative.
        ///
        ///Creates the derivative variable with the values of the function
        ///and first and second derivatives as specified. All higher order
        ///derivatives are zero.
        CubicSplineDerivatives(
            ///The value of the function (zeroth derivative)
            double func_value,

            ///The first derivative.
            double first_deriv,

            ///The second derivative.
            double second_deriv);

        ///Returns the derivative of the given order (zero is allowed).
        double order(unsigned deriv_order=1) const;
    };

    ///\brief A class representing a once differentiable function of a single
    ///argument
    class LIB_LOCAL OneArgumentDiffFunction
        : public OneArgumentFunction<double,double> {
    private :
#ifndef NO_SERIALIZE
        ///Needed for serialization to work.
        friend class boost::serialization::access;

        ///Serialize this function.
        ///The second parameter is supposed to be version
        template<class Archive>
            void serialize(Archive & ar, const unsigned int) {
                ar & boost::serialization::base_object<
                        OneArgumentFunction<double,double>
                     >(*this);
            }
#endif
    public:
        ///\brief Returns a pointer to the derivative of the function.
        ///
        ///The use of a pointer allows avoiding potentially expensive copy
        ///opertaions.
        virtual const FunctionDerivatives *deriv(double x) const=0;

        virtual ~OneArgumentDiffFunction() {}
    };

    ///\brief The derivatives of an identically zero quantity.
    ///
    ///\ingroup StellarSystem_group
    class LIB_LOCAL ZeroDerivatives : public FunctionDerivatives {
    private:
        double __zeroth_deriv;
    public:
        ///Create a derivative of an identically zero quantity.
        ZeroDerivatives(double zeroth_deriv=0) : __zeroth_deriv(zeroth_deriv)
        {}

        ///The deriv_order-th derivative.
        double order(unsigned order=1) const
        {return (order == 0 ? __zeroth_deriv : 0);}
    };

    ///A class representing a function that is identically zero.
    class LIB_LOCAL ZeroFunction : public OneArgumentDiffFunction {
    private :
#ifndef NO_SERIALIZE
        ///Needed for serialization to work.
        friend class boost::serialization::access;

        ///Serialize this function.
        ///The second parameter is supposed to be version
        template<class Archive>
            void serialize(Archive & ar, const unsigned int) {
                ar & boost::serialization::base_object< OneArgumentDiffFunction >(
                    *this
                );
            }
#endif

    public:
        ///See OneArgumentDiffFunction::deriv()
        virtual const FunctionDerivatives *deriv(double) const
        {return new ZeroDerivatives;}

        ///The function value (=0).
        double operator()(double) const {return 0;}

        ///The function is defined everywhere.
        double range_high() const {return Inf;}

        ///The function is defined everywhere
        double range_low() const {return -Inf;}

        ///\brief An iterator over the ages where the quantity takes the
        ///given y value.
        InterpSolutionIterator crossings(double =0) const;
    };

} //End Core namespace.

#endif

