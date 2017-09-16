/**\file
 * 
 * \brief Declares a dissipating body with two ConstPhaseLadDissipatingZone
 * components.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TWO_ZONE_BODY_H
#define __TWO_ZONE_BODY_H

#include "../ExponentialDecayDiffRotBody.h"
#include "../SaturatingSkumanichWindBody.h"

///A body with two zones.
class TwoZoneBody : public ExponentialDecayDiffRotBody,
					public SaturatingSkumanichWindBody {
private:
	///The top zone.
	DissipatingZone &__envelope,

					///The bottom zone.
					&__core;
public:
	///Create a body from two zones.
	TwoZoneBody(
			///The outer zone.
			DissipatingZone &envelope,

			///The inner zone.
			DissipatingZone &core,
			
			///The timescale over which the two zones are coupled.
			double coupling_timescale,

			///The strength of the angular momentum loss wind.
			double wind_strength,
			
			///The frequency at which the wind saturates
			double wind_saturation_frequency)
		: ExponentialDecayDiffRotBody(coupling_timescale),
		  SaturatingSkumanichWindBody(wind_strength,
				  					  wind_saturation_frequency),
		  __envelope(envelope), __core(core) {}

	unsigned number_zones() const {return 2;}

	DissipatingZone &zone(unsigned zone_index)
	{
		assert(zone_index<=1);
		return (zone_index==0 ? __envelope : __core);
	}

	const DissipatingZone &zone(unsigned zone_index) const
	{
		assert(zone_index<=1);
		return (zone_index==0 ? __envelope : __core);
	}

	virtual ~TwoZoneBody() {}
};

#endif
