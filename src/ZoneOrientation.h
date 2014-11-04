/**\file
 * 
 * \brief Declares a class for orientations of zones of DissipatingBody
 * objects.
 *
 * \ingroup StellarSystem_group
 */

/**\brief Orientations of zones of bodies in a binary system.
 *
 * \ingroup StellarSystem_group
 */

#include "DissipationQuantities.h"
#include "Eigen/Dense"
#include "Common.h"

class ZoneOrientation {
private:

	///The inclination of the zone relative to the orbit.
	double __inclination,

		   ///The periapsis of the orbit in the equatorial frame of the zone.
		   __periapsis;
public:
	ZoneOrientation(double inclination=NaN, double periapsis=NaN) :
		__inclination(inclination), __periapsis(periapsis) {}

	///Changes the zone orientation.
	void configure(double inclination, double periapsis)
	{__inclination=inclination; __periapsis=periapsis;}

	///The angle between the angular momenta of the zone and the orbit.
	double inclination() const {return __inclination;}

	///The argument of periapsis of this zone minus the reference zone's
	double periapsis() const {return __periapsis;}
};

///Transforms a vector betwen the coordinates systems of two zones.
Eigen::Vector3d zone_to_zone_transform(
		///The zone whose coordinate system the vectors are currently in.
		const ZoneOrientation &from_zone,

		///The zone whose coordinate system we want to transform the vectors
		///to.
		const ZoneOrientation &to_zone, 
		
		///The vector to transform.
		const Eigen::Vector3d &vector,

		///Derivatives with respect to inclination and periapsis can be
		///computed (assuming vector does not depend on these), in addition
		///to the regular transform. It is an error to request another
		///derivative.
		Dissipation::Derivative deriv=Dissipation::NO_DERIV,

		///If deriv is not NO_DERIV, derivatives can be computed with respect
		///to quantities of the from_zone (if this argument is true) or the
		///to_zone (if false).
		bool with_respect_to_from=false);

///Transforms the inclination and periapsis of a zone between references.
void transform_zone_orientation(
		///The zone whose orientation we wish to transform
		const ZoneOrientation &zone,

		///The reference frame in which we want to express the zone's
		///orientation in in the old reference frame.
		const ZoneOrientation &reference,

		///Overwritten by the inclination of zone in the new reference frame.
		double &inclination,

		///Overwritten by the periapsis of zone in the new reference frame.
		double &periapsis);


