#include "DiskPlanetSystem.h"

void DiskPlanetSystem::release_surface_spin()
{
	unsigned nzones=primary().number_zones();
	std::valarray<double> angmom(nzones), zeros(0.0, nzones-1);
	for(unsigned zone_ind=0; zone_ind<nzones; ++zone_ind)
		angmom[zone_ind]=primary().zone(zone_ind).angular_momentum();
	configure(age(), NaN, NaN, &(angmom[0]), &(zeros[0]), &(zeros[0]),
			  SINGLE);
}

void DiskPlanetSystem::add_secondary()
{
	unsigned nzones=primary().number_zones()+secondary().number_zones();
	std::valarray<double> angmom(nzones), inclination(nzones), 
		periapsis(nzones-1);
	unsigned zone_ind=0;
	for(short body_ind=0; body_ind<2; ++body_ind) {
		const DissipatingBody &body=(body_ind==0 ? primary() : secondary());
		for(unsigned body_zone_ind=0; body_zone_ind<body.number_zones();
				++body_zone_ind) {
			const DissipatingZone &zone=body.zone(body_zone_ind);
			angmom[zone_ind]=zone.angular_momentum();
			assert(zone.inclination()==0);
			assert(zone.periapsis()==0);
			inclination[zone_ind]=__initial_inclination;
			if(zone_ind) periapsis[zone_ind-1]=0;
			++zone_ind;
		}
	}
	configure(age(), __initial_semimajor, __initial_eccentricity,
			  &(angmom[0]), &(inclination[0]), &(periapsis[0]), BINARY);
}

DiskPlanetSystem::DiskPlanetSystem(DissipatingBody &body1,
		DissipatingBody &body2, double initial_semimajor,
		double initial_eccentricity, double initial_inclination,
		double disk_lock_frequency,  double disk_dissipation_age,
		double secondary_formation_age) :
	BinarySystem(body1, body2),
	__initial_semimajor(initial_semimajor),
	__initial_eccentricity(initial_eccentricity),
	__initial_inclination(initial_inclination),
	__disk_lock_frequency(disk_lock_frequency),
	__disk_dissipation_age(disk_dissipation_age),
	__secondary_formation_age(secondary_formation_age)
{
	if(initial_eccentricity<0 || initial_eccentricity>1) {
		std::ostringstream msg;
		msg << "Invalid initial eccentricity: " << initial_eccentricity
			<< " encountered in DiskPlanetSystem constructor!";
		throw Error::BadFunctionArguments(msg.str());
	}
	if(__initial_inclination<0 || __initial_inclination>M_PI) {
		std::ostringstream msg;
		msg << "Invalid initial inclination: " << initial_eccentricity
			<< " encountered in DiskPlanetSystem constructor, must be "
			   "between 0 and pi!";
		throw Error::BadFunctionArguments(msg.str());
	}
	if(__secondary_formation_age<__disk_dissipation_age) {
		std::ostringstream msg;
		msg << "Forming the secondary (at age= " << secondary_formation_age
			<< ") before the disk dissipates (at age="
			<< disk_dissipation_age << ") is not supported at this time.";
		throw Error::BadFunctionArguments(msg.str());
	}
	body1.set_surface_lock_frequency(__disk_lock_frequency);
}

void DiskPlanetSystem::reached_critical_age(double age)
{
	if(age==__disk_dissipation_age) release_surface_spin();
	if(age==__secondary_formation_age) add_secondary();
}

double DiskPlanetSystem::next_stop_age() const
{
	if(age()<__disk_dissipation_age) 
		return __disk_dissipation_age;
	else if(age()<__secondary_formation_age)
		return __secondary_formation_age;
	else return Inf;
}


