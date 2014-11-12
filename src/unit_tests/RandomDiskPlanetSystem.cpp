#include "RandomDiskPlanetSystem.h"

void RandomDiskPlanetSystem::create_system(EvolModeType evol_mode,
		bool lags_flip_sign)
{
	using namespace SystemParameters;
	std::valarray<double> angmom(__zones.size());
	for(unsigned i=0; i<__zones.size(); ++i) {
		__zones[i]=new ConstPhaseLagDissipatingZone(__lags[i],
				__parameters[FIRST_ZONE_INERTIA+i],
				__parameters[FIRST_ZONE_RADIUS+i],
				__parameters[FIRST_ZONE_MASS+i],
				lags_flip_sign,
				__parameters[FIRST_ZONE_INERTIA_DERIV+i],
				__parameters[FIRST_ZONE_RADIUS_DERIV+i],
				(i%2==0 ? 0 : __parameters[FIRST_CORE_MASS_DERIV+i/2]));
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
		bool match_primary_inclinations, bool zero_primary_inclinations,
		bool match_primary_periapses, bool match_secondary_inclinations, 
		bool zero_secondary_inclinations, bool match_secondary_periapses,
		bool zero_secondary_periapses, bool lags_flip_sign)
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
	__parameters[PRIMARY_RADIUS_DERIV]=
		uniform_rand(0, 100)*__parameters[PRIMARY_RADIUS];
	__parameters[SECONDARY_RADIUS]=
		uniform_rand(0, __parameters[SEMIMAJOR]/2.44
						*std::pow(__parameters[SECONDARY_MASS]
								  /__parameters[PRIMARY_MASS], 1.0/3.0));
	__parameters[SECONDARY_RADIUS_DERIV]=
		uniform_rand(0, 100)*__parameters[SECONDARY_RADIUS];
	__parameters[PRIMARY_CORE_RADIUS]=uniform_rand(0, 1)
									  *__parameters[PRIMARY_RADIUS];
	__parameters[PRIMARY_CORE_RADIUS_DERIV]=
		uniform_rand(0, 100)*__parameters[PRIMARY_CORE_RADIUS];
	__parameters[SECONDARY_CORE_RADIUS]=uniform_rand(0, 1)
									  *__parameters[SECONDARY_RADIUS];
	__parameters[SECONDARY_CORE_RADIUS_DERIV]=
		uniform_rand(0, 100)*__parameters[SECONDARY_CORE_RADIUS];
	__parameters[PRIMARY_CORE_MASS]=uniform_rand(0, 1)
									*__parameters[PRIMARY_MASS];
	__parameters[PRIMARY_CORE_MASS_DERIV]=
		uniform_rand(-100, 100)*__parameters[PRIMARY_CORE_MASS];
	__parameters[SECONDARY_CORE_MASS]=uniform_rand(0, 1)
									  *__parameters[SECONDARY_MASS];
	__parameters[SECONDARY_CORE_MASS_DERIV]=
		uniform_rand(-100, 100)*__parameters[SECONDARY_CORE_MASS];
	__parameters[PRIMARY_ENV_INERTIA]=
		std::pow(10.0, uniform_rand(-5, 0))
		*(__parameters[PRIMARY_MASS]-__parameters[PRIMARY_CORE_MASS])
		*std::pow(__parameters[PRIMARY_RADIUS], 2);
	__parameters[PRIMARY_ENV_INERTIA_DERIV]=
		uniform_rand(0, 100)*__parameters[PRIMARY_ENV_INERTIA];
	__parameters[PRIMARY_CORE_INERTIA]=
		std::pow(10.0, uniform_rand(-5, 0))
		*__parameters[PRIMARY_CORE_MASS]
		*std::pow(__parameters[PRIMARY_CORE_RADIUS], 2);
	__parameters[PRIMARY_CORE_INERTIA_DERIV]=
		uniform_rand(0, 100)*__parameters[PRIMARY_CORE_INERTIA];
	__parameters[SECONDARY_ENV_INERTIA]=
		std::pow(10.0, uniform_rand(-5, 0))
		*(__parameters[SECONDARY_MASS]-__parameters[SECONDARY_CORE_MASS])
		*std::pow(__parameters[SECONDARY_RADIUS], 2);
	__parameters[SECONDARY_ENV_INERTIA_DERIV]=
		uniform_rand(0, 100)*__parameters[SECONDARY_ENV_INERTIA];
	__parameters[SECONDARY_CORE_INERTIA]=
		std::pow(10.0, uniform_rand(-5, 0))
		*__parameters[SECONDARY_CORE_MASS]
		*std::pow(__parameters[SECONDARY_CORE_RADIUS], 2);
	__parameters[SECONDARY_CORE_INERTIA_DERIV]=
		uniform_rand(0, 100)*__parameters[SECONDARY_CORE_INERTIA];
	__parameters[PRIMARY_INCLINATION_ENV]=
		(evol_mode==BINARY && !zero_primary_inclinations
		 ? uniform_rand(0, M_PI)
		 : 0);
	__parameters[PRIMARY_INCLINATION_CORE]=
		(match_primary_inclinations || zero_primary_inclinations
		 ? __parameters[PRIMARY_INCLINATION_ENV]
		 : uniform_rand(0, M_PI));
	__parameters[SECONDARY_INCLINATION_ENV]=
		(zero_secondary_inclinations ? 0 : uniform_rand(0, M_PI));
	__parameters[SECONDARY_INCLINATION_CORE]=
		(match_secondary_inclinations || zero_secondary_inclinations
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
					(uniform_rand(0, 1)<0.2 ? 0 : 
					 std::pow(10.0, uniform_rand(-6, 0)));
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
	create_system(evol_mode, lags_flip_sign);
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

void add_zone_quantity(std::ostream &os, const std::string &name,
		const RandomDiskPlanetSystem &system, unsigned first_ind)
{
	os << "\t" << name << "=(";
	for(unsigned quant_ind=0; quant_ind<4; ++quant_ind) {
		os << system.quantity(
			static_cast<SystemParameters::Quantity>(first_ind+quant_ind));
		if(quant_ind==0 || quant_ind==2)
			os << ", ";
		else if(quant_ind==1) os << "), (";
	}
	os << ")" << std::endl;
}

void add_body_quantity(std::ostream &os, const std::string &name,
		const RandomDiskPlanetSystem &system, unsigned first_ind)
{
	os << "\t" << name << "=("
		<< system.quantity(
				static_cast<SystemParameters::Quantity>(first_ind))
		<< ", " << system.quantity(
				static_cast<SystemParameters::Quantity>(first_ind+1))
		<< ")" << std::endl;
}

std::ostream &operator<<(std::ostream &os,
		const RandomDiskPlanetSystem &system)
{
	using namespace SystemParameters;
	os << "\tt=" << system.quantity(AGE) << std::endl
	   << "\tWdisk=" << system.quantity(DISK_LOCK_FREQ) << std::endl
	   << "\tTdisk=" << system.quantity(DISK_DISSIPATION_AGE) << std::endl
	   << "\tTplanet=" << system.quantity(SECONDARY_FORMATION_AGE)
	   << std::endl
	   << "\ta=" << system.quantity(SEMIMAJOR) << std::endl
	   << "\te=" << system.quantity(ECCENTRICITY) << std::endl;
	add_zone_quantity(os, "zone masses", system, FIRST_ZONE_MASS);
	add_body_quantity(os, "core mass derivs", system, FIRST_CORE_MASS_DERIV);
	add_zone_quantity(os, "zone radii", system, FIRST_ZONE_RADIUS);
	add_zone_quantity(os, "zone radii derivs", system,
			FIRST_ZONE_RADIUS_DERIV);
	add_zone_quantity(os, "moments of inertia", system, FIRST_ZONE_INERTIA);
	add_zone_quantity(os, "moment of inertia derivs", system,
			FIRST_ZONE_INERTIA_DERIV);
	add_zone_quantity(os, "inclinations", system, FIRST_INCLINATION);
	add_zone_quantity(os, "periapses", system, FIRST_PERIAPSIS);
	add_zone_quantity(os, "angular velocities", system, FIRST_ZONE_ANGVEL);
	add_body_quantity(os, "wind strengths", system, FIRST_WIND_STRENGTH);
	add_body_quantity(os, "wind saturation freq.", system,
			FIRST_WIND_SAT_FREQ);
	add_body_quantity(os, "coupling timescales", system,
			FIRST_COUPLING_TIMESCALE);
	return os;
}
