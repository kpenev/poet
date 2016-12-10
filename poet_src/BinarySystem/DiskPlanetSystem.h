#ifndef __DISK_PLANET_SYSTEM_H
#define __DISK_PLANET_SYSTEM_H

#include "BinarySystem.h"
#include "Error.h"
#include <sstream>

/**\file
 *
 * \brief Declares a class of binary systems which start with a disk-locked
 * star which is then released and at some point in time a secondary appears.
 *
 * \ingroup StellarSystem_group
 */

/**\brief For some prescribed amount of time the surface of the pramary spins
 * at a prescribed rate, it is then released and (at a possibly different
 * age) a secondary body forms in a prescribed initial orbit.
 *
 * \ingroup StellarSystem_group
 */
class DiskPlanetSystem : virtual public BinarySystem {
private:
	///The semimajor axis of the orbit at which the secondary forms.
	double __initial_semimajor,

		   ///The eccentricity of the orbit at which the secondary forms.
		   __initial_eccentricity,

		   ///Inclination between surface zone of primary and initial orbit.
		   __initial_inclination,
		   
		   ///\brief Frequency of the surface spin of the primary when disk
		   ///is present.
		   __disk_lock_frequency,
		   
		   ///\brief Age when disk dissipates.
		   __disk_dissipation_age,
		   
		   ///\brief Age when the secondary forms.
		   __secondary_formation_age;

	///Releases the surface spin of the star when the disk dissipates.
	void release_surface_spin();

	///Adds the secondary to the system in its initial orbit.
	void add_secondary();
public:
	///Create the system.
	DiskPlanetSystem(
			///The first body in the system. Assumed to always be there, so
			///for a star-planet system this should be the star.
			DissipatingBody &body1, 

			///The second body in the system, must already have all its zones
			///configured as it will appear.
			DissipatingBody &body2,

			///The semimajor axis of the orbit at which the secondary forms.
			double initial_semimajor,

			///The eccentricity of the orbit at which the secondary forms.
			double initial_eccentricity,

			///Inclination between surface zone of primary and initial orbit.
			double initial_inclination,

			///\brief Frequency of the surface spin of the primary when disk
			///is present.
			double disk_lock_frequency,

			///\brief Age when disk dissipates.
			double disk_dissipation_age,

			///\brief Age when the secondary forms.
			double secondary_formation_age);

	///The age when the disk dissipates.
	double disk_dissipation_age() {return __disk_dissipation_age;}

	///\brief Change the system as necessary at the given age.
	///
	///Handles things like the disk dissipating and the secondary forming. 
	virtual void reached_critical_age(double age);

	///\brief The next age when the evolution needs to be stopped for a
	///system change
	virtual double next_stop_age() const;

	virtual ~DiskPlanetSystem() {}
};

#endif
