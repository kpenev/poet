#ifndef __SATURATING_SKUMANICH_WIND_BODY_H
#define __SATURATING_SKUMANICH_WIND_BODY_H

#include "DissipatingBody.h"
#include "BinarySystem.h"

/**\file
 *
 * \brief Decrales a body subject to angular momentum loss 
 * \f$\propto\omega\min(\omega, \omega_{sat})^2\f$.
 *
 * \ingroup StellarSystem_group
 */

///\brief A DissipatingBody which loses angular momentum at a rate
/// \f$\propto\omega\min(\omega, \omega_{sat})^2\f$
///
///\ingroup StellarSystem_group
class SaturatingSkumanichWindBody : virtual public DissipatingBody {
private:
	///The strength of the magnetic wind
	double __wind_strength,

		///The frequency at which the wind loss saturates in rad/day.
		__saturation_freq;

	///Is the wind currently saturated?
	bool __saturated;

	///The saturation states recorded by add_to_evolution() so far.
	std::list<bool> __saturation_evolution;
public:
	SaturatingSkumanichWindBody(
			///The strength of the wind.
			double wind_strength,

			///The frequency at which the wind loss saturates in rad/day.
			double saturation_frequency)
		: __wind_strength(wind_strength), 
		  __saturation_freq(saturation_frequency) {}

	///See DissipatingBody::angular_momentum_loss().
	double angular_momentum_loss(
			Dissipation::Derivative deriv=Dissipation::NO_DERIV) const;

	///The frequency at which the wind loss saturates in rad/day.
	double saturation_frequency() {return __saturation_freq;}

	///Sets the saturation based on the currently configured spin frequency.
	void detect_saturation() 
	{__saturated=(spin_frequency()>__saturation_freq);}

	///Is the wind loss currently saturated?
	bool saturated() {return __saturated;}

	///Called by the stopping condition monitoring wind saturation.
	void saturation_freq_crossed(
			///The sign of the rate of change of the spin frequency when it
			///was equal to the saturation frequency.
			short
#ifdef DEBUG
			deriv_sign
#endif
			)
	{
#ifdef DEBUG
		assert(deriv_sign==(__saturated ? -1 : 1));
#endif
		__saturated=!__saturated;
	}

	///Appends the state defined by last configure(), to the evolution.
	virtual void add_to_evolution()
	{
		__saturation_evolution.push_back(__saturated); 
		DissipatingBody::add_to_evolution();
	}

	///Discards the last steps from the evolution.
	virtual void rewind_evolution(
			///How many steps of evolution to discard.
			unsigned nsteps)
	{
		for(unsigned i=0; i<nsteps; ++i) __saturation_evolution.pop_back();
		DissipatingBody::rewind_evolution(nsteps);
	}

	///Discards all evolution
	virtual void reset_evolution()
	{
		__saturation_evolution.clear();
		DissipatingBody::reset_evolution();
	}

	///\brief Conditions detecting the next possible discontinuities in the
	///evolution due to this body.
	///
	///Must be deleted when no longer necessary.
	virtual CombinedStoppingCondition *stopping_conditions(
			///The system being evolved.
			BinarySystem &system, 

			///Is the body the primary in the system.
			bool primary);

	///The tabulated wind saturation states so far.
	const std::list<bool> &wind_saturation_evolution() const
	{return __saturation_evolution;}

	///Resets its saturation state after a discontinous spin jump.
	void spin_jumped()
	{
		detect_saturation();
		DissipatingBody::spin_jumped();
	}
};

#endif
