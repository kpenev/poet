#ifndef __RANDOM_DISK_PLANET_SYSTEM_H
#define __RANDOM_DISK_PLANET_SYSTEM_H

#include "../DiskPlanetSystem.h"
#include "TwoZoneBody.h"
#include "ConstPhaseLagDissipatingZone.h"
#include "Common.h"

namespace SystemParameters {
	enum Quantity {
		AGE,
		DISK_LOCK_FREQ,
		DISK_DISSIPATION_AGE,
		SECONDARY_FORMATION_AGE,
		SEMIMAJOR,
		ECCENTRICITY,
		FIRST_ZONE_MASS,
		PRIMARY_MASS=FIRST_ZONE_MASS,
		PRIMARY_CORE_MASS,
		SECONDARY_MASS,
		SECONDARY_CORE_MASS,
		FIRST_CORE_MASS_DERIV,
		PRIMARY_CORE_MASS_DERIV=FIRST_CORE_MASS_DERIV,
		SECONDARY_CORE_MASS_DERIV,
		FIRST_ZONE_RADIUS,
		PRIMARY_RADIUS=FIRST_ZONE_RADIUS,
		PRIMARY_CORE_RADIUS,
		SECONDARY_RADIUS,
		SECONDARY_CORE_RADIUS,
		FIRST_ZONE_RADIUS_DERIV,
		PRIMARY_RADIUS_DERIV=FIRST_ZONE_RADIUS_DERIV,
		PRIMARY_CORE_RADIUS_DERIV,
		SECONDARY_RADIUS_DERIV,
		SECONDARY_CORE_RADIUS_DERIV,
		FIRST_ZONE_INERTIA,
		PRIMARY_ENV_INERTIA=FIRST_ZONE_INERTIA,
		PRIMARY_CORE_INERTIA,
		SECONDARY_ENV_INERTIA,
		SECONDARY_CORE_INERTIA,
		FIRST_ZONE_INERTIA_DERIV,
		PRIMARY_ENV_INERTIA_DERIV=FIRST_ZONE_INERTIA_DERIV,
		PRIMARY_CORE_INERTIA_DERIV,
		SECONDARY_ENV_INERTIA_DERIV,
		SECONDARY_CORE_INERTIA_DERIV,
		FIRST_INCLINATION,
		PRIMARY_INCLINATION_ENV=FIRST_INCLINATION,
		PRIMARY_INCLINATION_CORE,
		SECONDARY_INCLINATION_ENV,
		SECONDARY_INCLINATION_CORE,
		FIRST_PERIAPSIS,
		PRIMARY_PERIAPSIS_CORE=FIRST_PERIAPSIS,
		SECONDARY_PERIAPSIS_ENV,
		SECONDARY_PERIAPSIS_CORE,
		FIRST_ZONE_ANGVEL,
		PRIMARY_ANGVEL_ENV=FIRST_ZONE_ANGVEL,
		PRIMARY_ANGVEL_CORE,
		SECONDARY_ANGVEL_ENV,
		SECONDARY_ANGVEL_CORE,
		FIRST_WIND_STRENGTH,
		PRIMARY_WIND_STRENGTH=FIRST_WIND_STRENGTH,
		SECONDARY_WIND_STRENGTH,
		FIRST_WIND_SAT_FREQ,
		PRIMARY_WIND_SAT_FREQ=FIRST_WIND_SAT_FREQ,
		SECONDARY_WIND_SAT_FREQ,
		FIRST_COUPLING_TIMESCALE,
		PRIMARY_COUPLING_TIMESCALE=FIRST_COUPLING_TIMESCALE,
		SECONDARY_COUPLING_TIMESCALE,
		NUM_QUANTITIES
	};
};

class RandomDiskPlanetSystem {
private:
	///The values of the randomly selected parameters.
	std::valarray<double> __parameters;

	///The lags for each zone.
	std::vector<Lags> __lags;

	///The spin orbit lock for each zone (not all will be locked in general).
	std::vector<SpinOrbitLockInfo> __locks;

	///How many zones were actually locked.
	unsigned __num_locked_zones;

	///The four zones in the body.
	std::vector<ConstPhaseLagDissipatingZone *> __zones;

	///The primary and secondary body in the system.
	std::vector<TwoZoneBody *> __bodies;

	///The system.
	DiskPlanetSystem *__system;

	///The orbital frequnecy.
	double __worb;

	///Uses the contents of __parameters to configure the system.
	void create_system(EvolModeType evol_mode, bool lags_flip_sign);

	///Locks a random number of zones in random ratios with the orbit.
	void lock_zones(
			///The minimum number of zones to lock.
			unsigned min_locked_zones,

			///The maximum number of zones to lock.
			unsigned max_locked_zones);
public:
	///Generates and configures a system randomly.
	RandomDiskPlanetSystem(EvolModeType evol_mode,
			unsigned min_locked_zones=0,
			unsigned max_locked_zones=4,
			bool circular=false,
			bool match_primary_inclinations=false,
			bool zero_primary_inclinations=false,
			bool match_primary_periapses=false,
			bool match_secondary_inclinations=false,
			bool zero_secondary_inclinations=false,
			bool match_secondary_periapses=false,
			bool zero_secondary_periapses=false,
			bool lags_flip_sign=false);

	///Returns the random value chosen for a quantity.
	double quantity(SystemParameters::Quantity q) const
	{return __parameters[q];}

	///The lag assigned to a zone.
	const Lags &lags(
			///Which zone's lag? Indices go from primary envelope to
			///secondary core.
			unsigned zone_ind) const {return __lags[zone_ind];}

	///Returns a reference to a locally held system.
	DiskPlanetSystem &operator()() {return *__system;}

	///Returns a reference to a locally held system.
	const DiskPlanetSystem &operator()() const {return *__system;}

	///The lock information of a single zone.
	const SpinOrbitLockInfo &lock(unsigned zone_ind) const
	{return __locks[zone_ind];}

	///The number of locked zones.
	unsigned num_locked_zones() const {return __num_locked_zones;}

	///Returns the orbital frequency of the generated system.
	double orbital_frequency() const {return __worb;}

	///Cleans up.
	~RandomDiskPlanetSystem();
};

///Describe the system on the given stream.
std::ostream &operator<<(std::ostream &os,
		const RandomDiskPlanetSystem &system);

#endif
