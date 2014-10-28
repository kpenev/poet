#ifndef __SATURATING_SKUMANICH_WIND_BODY_H
#define __SATURATING_SKUMANICH_WIND_BODY_H

#include "DissipatingBody.h"

/**\file
 *
 * \brief Decrales a body subject to angular momentum loss 
 * \f$\propt\omega\min(\omega, \omega_{sat})^2\f$.
 *
 * \ingroup StellarSystem_group
 */

///\brief A DissipatingBody which loses angular momentum at a rate
/// \f$\propt\omega\min(\omega, \omega_{sat})^2\f$
///
///\ingroup StellarSystem_group
class SaturatingSkumanichWindBody : public DissipatingBody {
private:
	///The frequency at which the wind loss saturates in rad/day.
	double __saturation_freq;

	///Is the wind currently saturated?
	bool __saturated;
public:
	SaturatingSkumanichWindBody(
			///The frequency at which the wind loss saturates in rad/day.
			double saturation_frequency)
		: __saturation_freq(saturation_frequency) {}

	///The frequency at which the wind loss saturates in rad/day.
	double saturation_frequency() {return __saturation_freq;}

	///Is the wind loss currently saturated?
	bool saturated() {return __saturated;}

	///Called by the stopping condition monitoring wind saturation.
	void saturation_freq_crossed(
			///The sign of the rate of change of the spin frequency when it
			///was equal to the saturation frequency.
			short deriv_sign)
	{
#ifdef DEBUG
		assert(deriv_sign==(__saturated ? -1 : 1));
#endif
		__saturated=!__saturated;
	}
};

#endif
