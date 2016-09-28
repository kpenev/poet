/**\file
 *
 * \brief Defines the class that calculates \f$I_{2m,2n}\f$ from the 
 * description.
 *
 * \ingroup Utilities
 */

#ifndef __SIN_COS_INTEGRALS_H
#define __SIN_COS_INTEGRALS_H

#include <vector>
#include <cmath>

///\brief Calculates \f$\int_0^{2\pi} \sin^{2m} u \cos^{2n} u du\f$.
///
///\ingroup Utilities
class SinCosIntegrals {
private:
	///Amount of space to reserve for newly created m entries in __values.
	unsigned __pre_allocate;
	
	///\brief The currently computed values.
	///
	///The inner index is n and the outer in m.
	std::vector< std::vector<double> > __values;

	///
public:
	///Initialize an object with only the m=n=0 term calculated.
	SinCosIntegrals(unsigned pre_allocate=0);

	///Returns the value of \f$\int_0^{2\pi} \sin^{2m} u \cos^{2n} u du\f$.
	double operator()(unsigned m, unsigned n);
};

#endif
