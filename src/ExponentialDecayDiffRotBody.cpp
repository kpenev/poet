#include "ExponentialDecayDiffRotBody.h"

Eigen::Vector3d ExponentialDecayDiffRotBody::angular_momentum_coupling(
		unsigned top_zone_index, Dissipation::Derivative deriv,
		bool with_respect_to_top) const
{
#ifdef DEBUG
	assert(top_zone_index<number_zones()-1);
#endif
	if(		deriv==Dissipation::AGE ||
			deriv==Dissipation::ORBITAL_FREQUENCY ||
			deriv==Dissipation::ECCENTRICITY ||
			deriv==Dissipation::RADIUS ||
			deriv==Dissipation::SEMIMAJOR)
		return Eigen::Vector3d(0, 0, 0);
	DissipatingZone &top_zone=zone(top_zone_index),
					&bottom_zone=zone(top_zone_index+1);
	double top_inertia=top_zone.moment_of_inertia(),
		   bot_inertia=bottom_zone.moment_of_inertia(),
		   sum_inertia=top_inertia+bot_inertia;
	else if(deriv==Dissipation::NO_DERIV || 
			deriv==Dissipation::MOMENT_OF_INERTIA) {
		Eigen::Vector3d bottom_angmom=zone_to_zone_transform(
				bottom_zone, top_zone, 
				Eigen::Vector3d(0, 0, bottom_zone.angular_momentum())),
				angmom_diff=(top_inertia*bottom_angmom
							 -
							 Eigen::Vector3d(0, 0, 
								 bot_inertia()*top_zone.angular_momentum()))
				   			/sum_inertia;
		if(deriv==Dissipation::NO_DERIV) 
			return angmom_diff/__coupling_timescale;
		else if(with_respect_to_top) 
			return (bottom_angmom-angmom_diff)/__coupling_timescale;
		else {
			angmom_diff[2]+=top_zone.angular_momentum();
			return (-1.0/__coupling_timescale)*angmom_diff;
		}
	} else if(deriv==Dissipation::INCLINATION || 
			  deriv==Dissipation::PERIAPSIS) 
		return top_inertia()/sum_inertia/__coupling_timescale
			   *
			   zone_to_zone_transform(bottom_zone, top_zone, 
					   Eigen::Vector3d(0, 0, bottom_zone.angular_momentum()),
					   deriv, with_respect_to_top);
	else if(deriv==Dissipation::SPIN_FREQUENCY) {
		double scale=
			top_inertia*bot_inertia/sum_inertia/__coupling_timescale;
		if(with_respect_to_top) return Eigen::Vector3d(0, 0, scale);
		else return zone_to_zone_transform(bottom_zone, top_zone,
				Eigen::Vector3d(0, 0, scale));
	} else {
#ifdef DEBUG
		assert(deriv==Dissipation::SPIN_ANGMOM);
#endif
		if(with_respect_to_top) 
			return Eigen::Vector3d(0, 0, 
					-bot_inertia/sum_inertia/__coupling_timescale);
		else return zone_to_zone_transform(bottom_zone, top_zone,
				Eigen::Vector3d(0, 0, top_inertia/sum_inertia))
			/__coupling_timescale;
	}
}
