#include "RandomDiskPlanetSystem.h"

void RandomDiskPlanetSystem::create_system(EvolModeType evol_mode)
{
	using namespace SystemParameters;
	std::valarray<double> angmom(__zones.size());
	for(unsigned i=0; i<__zones.size(); ++i) {
		__zones[i]=new ConstPhaseLagDissipatingZone(__lags[i],
				__parameters[FIRST_ZONE_INERTIA+i],
				__parameters[FIRST_ZONE_RADIUS+i],
				__parameters[FIRST_ZONE_MASS+i]);
		angmom[i]=__parameters[FIRST_ZONE_INERTIA+i]
				  *__parameters[FIRST_ZONE_ANGVEL+i];
	}
	for(unsigned i=0; i<__bodies.size(); ++i)
		__bodies[i]=new TwoZoneBody(*__zones[2*i], *__zones[2*i+1],
									__parameters[FIRST_COUPLING_TIMESCALE+i],
									__parameters[FIRST_WIND_STRENGTH+i],
									__parameters[FIRST_WIND_SAT_FREQ+i]);
	__system=new DiskPlanetSystem(*__bodies[0], *__bodies[1],
								  __parameters[SEMIMAJOR],
								  __parameters[ECCENTRICITY],
								  __parameters[FIRST_INCLINATION],
								  __parameters[DISK_LOCK_FREQ],
								  __parameters[DISK_DISSIPATION_AGE],
								  __parameters[SECONDARY_FORMATION_AGE]);
	__system->configure(__parameters[AGE],
					    (evol_mode==BINARY ? __parameters[SEMIMAJOR] : NaN),
					    (evol_mode==BINARY ? __parameters[ECCENTRICITY]
						 				   : NaN),
					    &(angmom[(evol_mode==LOCKED_SURFACE_SPIN ? 1 : 0)]),
					    (evol_mode==LOCKED_SURFACE_SPIN
						 ? NULL
						 : &(__parameters[(evol_mode==SINGLE ? 1 : 0)
							 			  +FIRST_INCLINATION])),
					    (evol_mode==LOCKED_SURFACE_SPIN 
						 ? NULL
						 : &(__parameters[FIRST_PERIAPSIS])),
					    evol_mode);
	__bodies[0]->detect_saturation();
	__bodies[1]->detect_saturation();
}

void RandomDiskPlanetSystem::lock_zones(unsigned min_locked_zones,
		unsigned max_locked_zones)
{
	__num_locked_zones=std::rand()%(max_locked_zones + 1 - min_locked_zones)
					   +
					   min_locked_zones;
	std::list<unsigned> still_unlocked;
	for(unsigned i=0; i<4; ++i) still_unlocked.push_back(i);
	for(unsigned i=0; i<__num_locked_zones; ++i) {
		std::list<unsigned>::iterator unlocked_i=still_unlocked.begin();
		for(unsigned to_lock=std::rand()%still_unlocked.size(); to_lock>0;
				--to_lock) ++unlocked_i;
		int orb_freq_mult=rand()%5-3, spin_freq_mult=rand()%2+1;
		__bodies[*unlocked_i/2]->lock_zone_spin(*unlocked_i%2, orb_freq_mult,
												spin_freq_mult);
		__parameters[SystemParameters::FIRST_ZONE_ANGVEL+*unlocked_i]=
			__worb*static_cast<double>(orb_freq_mult)
			/static_cast<double>(spin_freq_mult);
		__locks[*unlocked_i].set_lock(orb_freq_mult, spin_freq_mult);
		still_unlocked.erase(unlocked_i);
	}
}

RandomDiskPlanetSystem::RandomDiskPlanetSystem(EvolModeType evol_mode,
		unsigned min_locked_zones, unsigned max_locked_zones, bool circular,
		bool match_primary_inclinations, bool match_primary_periapses, 
		bool match_secondary_inclinations, bool match_secondary_periapses,
		bool zero_secondary_periapses)
	: __parameters(SystemParameters::NUM_QUANTITIES), __lags(4),
	__locks(4), __zones(4, NULL), __bodies(2, NULL), __system(NULL)
{
	using namespace SystemParameters;
	__parameters[AGE]=uniform_rand(0, 10);
	__parameters[SEMIMAJOR]=uniform_rand(3, 30);
	__parameters[ECCENTRICITY]=(circular ? 0 : uniform_rand(0, 1));
	__parameters[PRIMARY_MASS]=std::pow(10.0, uniform_rand(-3.0, 2.0));
	__parameters[SECONDARY_MASS]=std::pow(10.0, uniform_rand(-3.0, 0.0))
								 *__parameters[PRIMARY_MASS];
	__parameters[PRIMARY_RADIUS]=uniform_rand(0, 1)*__parameters[SEMIMAJOR];
	__parameters[SECONDARY_RADIUS]=
		uniform_rand(0, __parameters[SEMIMAJOR]/2.44
						*std::pow(__parameters[SECONDARY_MASS]
								  /__parameters[PRIMARY_MASS], 1.0/3.0));
	__parameters[PRIMARY_CORE_RADIUS]=uniform_rand(0, 1)
									  *__parameters[PRIMARY_RADIUS];
	__parameters[SECONDARY_CORE_RADIUS]=uniform_rand(0, 1)
									  *__parameters[SECONDARY_RADIUS];
	__parameters[PRIMARY_CORE_MASS]=uniform_rand(0, 1)
									*__parameters[PRIMARY_MASS];
	__parameters[SECONDARY_CORE_MASS]=uniform_rand(0, 1)
									  *__parameters[SECONDARY_MASS];
	__parameters[PRIMARY_ENV_INERTIA]=
		std::pow(10.0, uniform_rand(-5, 0))
		*(__parameters[PRIMARY_MASS]-__parameters[PRIMARY_CORE_MASS])
		*std::pow(__parameters[PRIMARY_RADIUS], 2);
	__parameters[PRIMARY_CORE_INERTIA]=
		std::pow(10.0, uniform_rand(-5, 0))
		*__parameters[PRIMARY_CORE_MASS]
		*std::pow(__parameters[PRIMARY_CORE_RADIUS], 2);
	__parameters[SECONDARY_ENV_INERTIA]=
		std::pow(10.0, uniform_rand(-5, 0))
		*(__parameters[SECONDARY_MASS]-__parameters[SECONDARY_CORE_MASS])
		*std::pow(__parameters[SECONDARY_RADIUS], 2);
	__parameters[SECONDARY_CORE_INERTIA]=
		std::pow(10.0, uniform_rand(-5, 0))
		*__parameters[SECONDARY_CORE_MASS]
		*std::pow(__parameters[SECONDARY_CORE_RADIUS], 2);
	__parameters[PRIMARY_INCLINATION_ENV]=
		(evol_mode==BINARY ? uniform_rand(0, M_PI) : 0);
	__parameters[PRIMARY_INCLINATION_CORE]=
		(match_primary_inclinations ? __parameters[PRIMARY_INCLINATION_ENV]
		 							: uniform_rand(0, M_PI));
	__parameters[SECONDARY_INCLINATION_ENV]=uniform_rand(0, M_PI);
	__parameters[SECONDARY_INCLINATION_CORE]=
		(match_secondary_inclinations
		 ? __parameters[SECONDARY_INCLINATION_ENV]
		 : uniform_rand(0, M_PI));
	__parameters[PRIMARY_PERIAPSIS_CORE]=(match_primary_periapses
										  ? 0
										  : uniform_rand(0, 2.0*M_PI));
	__parameters[SECONDARY_PERIAPSIS_ENV]=(zero_secondary_periapses
										   ? 0
										   : uniform_rand(0, 2.0*M_PI));
	__parameters[SECONDARY_PERIAPSIS_CORE]=
		(zero_secondary_periapses || match_secondary_periapses
		 ? __parameters[SECONDARY_PERIAPSIS_ENV]
		 : uniform_rand(0, 2.0*M_PI));
	double __worb=orbital_angular_velocity(__parameters[PRIMARY_MASS],
										   __parameters[SECONDARY_MASS],
										   __parameters[SEMIMAJOR]);
	__parameters[DISK_LOCK_FREQ]=std::pow(10.0, uniform_rand(-1, 1))*__worb;
	switch(evol_mode) {
		case LOCKED_SURFACE_SPIN :
			__parameters[DISK_DISSIPATION_AGE]=
				uniform_rand(__parameters[AGE], 10);
			__parameters[SECONDARY_FORMATION_AGE]=
				uniform_rand(__parameters[DISK_DISSIPATION_AGE], 10);
			break;
		case BINARY :
			__parameters[DISK_DISSIPATION_AGE]=
				uniform_rand(0, __parameters[AGE]);
			__parameters[SECONDARY_FORMATION_AGE]=
				uniform_rand(__parameters[DISK_DISSIPATION_AGE],
							 __parameters[AGE]);
			break;
		case SINGLE :
			__parameters[DISK_DISSIPATION_AGE]=
				uniform_rand(0, __parameters[AGE]);
			__parameters[SECONDARY_FORMATION_AGE]=
				uniform_rand(__parameters[AGE], 10);
			break;
		default : assert(false);
	}
	for(unsigned i=0; i<__zones.size(); ++i) {
		__parameters[FIRST_ZONE_ANGVEL+i]=
			std::pow(10.0, uniform_rand(-1, 1))*__worb;
		for(int m=-2; m<=2; ++m)
			for(int mp=-2; mp<=0; ++mp) {
				__lags[i](m, mp)=
					(uniform_rand(0, 1)<0.2 ? 0 : uniform_rand(0, 10));
				__lags[i](-m, -mp)=-__lags[i](m, mp);
			}
	}
	for(unsigned i=0; i<__bodies.size(); ++i) {
		__parameters[FIRST_WIND_STRENGTH+i]=uniform_rand(0,1);
		__parameters[FIRST_WIND_SAT_FREQ+i]=
			std::pow(10.0, uniform_rand(-1, 1))*__worb;
		__parameters[FIRST_COUPLING_TIMESCALE+i]=
			std::pow(10.0, uniform_rand(-6,1));
	}
	create_system(evol_mode);
	if(evol_mode==BINARY) lock_zones(min_locked_zones, max_locked_zones);
}

RandomDiskPlanetSystem::~RandomDiskPlanetSystem()
{
	for(unsigned i=0; i<__zones.size(); ++i)
		if(__zones[i]) delete __zones[i];
	for(unsigned i=0; i<__bodies.size(); ++i)
		if(__bodies[i]) delete __bodies[i];
	if(__system) delete __system;
}


