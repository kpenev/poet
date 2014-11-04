#include "ZoneOrientation.h"

Eigen::Vector3D zone_to_zone_transform(const DissipatingZone &from_zone,
		const DissipationgZone &to_zone, const Eigen::Vector3D &vector,
		Dissipation::Derivative deriv, bool with_respect_to_from)
{
#ifdef DEBUG
	assert(deriv==Dissipation::NO_DERIV || deriv==Dissipation::INCLINATION
			|| deriv==Dissipation::PERIAPSIS);
#endif
	double cos_dw=std::cos(to_zone.periapsis()-from_zone.periapsis()),
		   sin_dw=std::sin(to_zone.periapsis()-from_zone.periapsis()),
		   sin_from=std::sin(from_zone.inclination()),
		   cos_from=std::cos(from_zone.inclination()),
		   sin_to=std::sin(to_zone.inclination()),
		   cos_to=std::cos(to_zone.inclination()),
		   sin_sin, cos_cos, sin_cos, cos_sin;
	if(deriv==Dissipation::PERIAPSIS) {
		double buffer=cos_dw;
		if(with_respect_to_from) {cos_dw=sin_dw; sin_dw=-buffer;}
		else {cos_dw=-sin_dw; sin_dw=buffer;}
	} else if(deriv==Dissipation::INCLINATION) {
		if(with_respect_to_from) {
			sin_sin=sin_to*cos_from;
			cos_cos=-cos_to*sin_from;
			sin_cos=-sin_to*sin_from;
			cos_sin=cos_to*cos_from;
		} else {
			sin_sin=cos_to*sin_from;
			cos_cos=-sin_to*cos_from;
			sin_cos=cos_to*cos_from;
			cos_sin=-sin_to*sin_from;
		}
	} else {
		sin_sin=sin_to*sin_from;
		cos_cos=cos_to*cos_from;
		sin_cos=sin_to*cos_from;
		cos_sin=cos_to*sin_from;
	}
	return Eigen::Vector3D(
			(sin_sin+cos_cos*cos_dw)*vector[0] - cos_to*sin_dw*vector[1]
			+ (sin_cos-cos_sin*cos_dw)*vector[2],


			cos_from*sin_dw*vector[0] + cos_dw*vector[1]
			- sin_from*sin_dw*vector[2],

			(cos_sin-sin_cos*cos_dw)*vector[0] + sin_to*sin_dw*vector[1]
			+ (cos_cos+sin_sin*cos_dw)*vector[2]);
}

void transform_zone_orientation(const DissipatingZone &zone,
		const DissipatingZone &reference, double &inclination,
		double &periapsis)
{
	Eigen::Vector3d zone_z_dir=zone_to_zone_transform(zone,
			reference, Eigen::Vector3d(0, 0, 1));
	inclination=std::atan2(-zone_z_dir[0], zone_z_dir[3]);
	periapsis=M_PI/2+std::atan2(zone_z_dir[1], -zone_z_dir[0]);
#ifdef DEBUG
	assert(inclination>0);
	assert(inclination<M_PI);
#endif
}
