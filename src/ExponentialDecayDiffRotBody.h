#ifndef __EXPONENTIAL_DECAY_DIFF_ROT_BODY_H
#define __EXPONENTIAL_DECAY_DIFF_ROT_BODY_H

#include "DissipatingBody.h"

class ExponentialDecayDiffRotBody : public DissipatingBody {
private:
	///The 
	double __coupling_timescale;
public:
	///\brief Create a DissipatingBody whose differential rotation between
	///decays exponentially.
	ExponentialDecayDiffRotBody(
			///The e-folding timsecale for the differential rotation decay.
			double coupling_timescale)
		: __coupling_timescale(coupling_timescale) {}

	///\brief Coupling torque for two neighboring zones in the coordinate
	///system of the top zone.
	///
	///Units:
	/// \f$M_\odot R_\odot^2\mathrm{rad}\mathrm{day}^{-1}mathrm{Gyr}^{-1}\f$.
	Eigen::Vector3d angular_momentum_coupling(
			///The index of the outer of the two zones whose coupling is
			///needed.
			unsigned top_zone_index,
			
			///The derivative of the coupling torque required.
			Dissipation::Derivative deriv=Dissipation::NO_DERIV,
			
			///For derivatives with respect to zone specific quantities,
			///this determines which zone's quantity to differentiate with
			///respect to (top zone if true, bottom zone if false).
			bool with_respect_to_top=false) const;
}

#endif
