#include "ConstSolutionIterator.h"

std::list<double> ConstSolutionIterator::__placeholder_list;
std::list<double>::const_iterator 
	ConstSolutionIterator::__placeholder_iterator=__placeholder_list.end();

void ConstSolutionIterator::create_missing_lists(
		const std::list<double> &tabulation_ages
)
{
	const DissipatingZone &convective=__star.zone(0),
		  				  &radiative=__star.zone(1);
	static std::vector< std::list<double> >
		quantities(OutCol::LAST_NO_ORBIT+1);
	for(
			std::list<double>::const_iterator t_i=tabulation_ages.begin();
			t_i!=tabulation_ages.end();
			++t_i
	) {
		double age=*t_i;
		quantities[OutCol::ICONV].push_back(
				convective.moment_of_inertia(age)
		);
		quantities[OutCol::IRAD].push_back(
				radiative.moment_of_inertia(age)
		);
		quantities[OutCol::RSTAR].push_back(convective.outer_radius(age));
		quantities[OutCol::RRAD].push_back(radiative.outer_radius(age));
		quantities[OutCol::MRAD].push_back(radiative.outer_mass(age));
		quantities[OutCol::ICONV_DERIV].push_back(
				convective.moment_of_inertia(age, 1)
		);
		quantities[OutCol::IRAD_DERIV].push_back(
				radiative.moment_of_inertia(age, 1)
		);
		quantities[OutCol::RSTAR_DERIV].push_back(
				convective.outer_radius(age, 1)
		);
		quantities[OutCol::RRAD_DERIV].push_back(
				radiative.outer_radius(age, 1)
		);
		quantities[OutCol::MRAD_DERIV].push_back(
				radiative.outer_mass(age, 1)
		);
		quantities[OutCol::ICONV_SECOND_DERIV].push_back(
				convective.moment_of_inertia(age, 2)
		);
		quantities[OutCol::IRAD_SECOND_DERIV].push_back(
				radiative.moment_of_inertia(age, 2)
		);
		quantities[OutCol::RRAD_SECOND_DERIV].push_back(
				radiative.outer_radius(age, 2)
		);
	}
	for(unsigned i=0; i<=OutCol::LAST_NO_ORBIT; ++i)
		if(i!=OutCol::AGE) __real_iterators[i]=quantities[i].begin();
}

void ConstSolutionIterator::fix_no_evolution(
		double start_age,
		double end_age,
		double timestep,
		const std::list<double> &required_ages
)
{
	static std::list<double> age_list;
	std::list<double>::const_iterator
		required_ages_iter=required_ages.begin();
	for(double age=start_age; age<=end_age; age+=timestep) {
		for(;required_ages_iter!=required_ages.end() &&
				*required_ages_iter<age; ++required_ages_iter)
			age_list.push_back(*required_ages_iter);
		age_list.push_back(age);
	}
	__real_iterators[OutCol::AGE]=age_list.begin();
	__last_age=age_list.end();
	create_missing_lists(age_list);
	for(
			unsigned i=OutCol::LAST_NO_ORBIT+1;
			i<OutCol::NUM_REAL_OUTPUT_QUANTITIES;
			++i
	) __real_iterators[i]=__placeholder_iterator;
	for(
			unsigned i=0;
			i!=OutCol::NUM_REAL_OUTPUT_QUANTITIES;
			++i
	) std::cerr << output_column_names()[i] << " is valid: "
				<< (__real_iterators[i]!=__placeholder_iterator)
				<< std::endl;
}

ConstSolutionIterator::ConstSolutionIterator(
		const OrbitSolver &solver,
		const BinarySystem &system,
		const YRECStar &star,
		double start_age,
		double end_age,
		double timestep,
		const std::list<double> &required_ages
) : __real_iterators(OutCol::NUM_REAL_OUTPUT_QUANTITIES,
					 __placeholder_iterator),
	__mode(solver.mode_evolution().begin()),
	__wind_saturation(star.wind_saturation_evolution().begin()),
	__last_age(solver.evolution_ages().end()),
	__mstar(star.mass()),
	__mplanet(system.secondary().mass()),
	__star(star)
{
	const DissipatingZone &convective=star.zone(0),
		  				  &radiative=star.zone(1);
	__real_iterators[OutCol::AGE]=solver.evolution_ages().begin();
	__real_iterators[OutCol::ICONV]=convective.get_evolution_real(
			MOMENT_OF_INERTIA_FIRST_DERIV
	).begin();
	__real_iterators[OutCol::IRAD]=radiative.get_evolution_real(
			MOMENT_OF_INERTIA_FIRST_DERIV
	).begin();
	__real_iterators[OutCol::RSTAR]=convective.get_evolution_real(
			OUTER_RADIUS
	).begin();
	__real_iterators[OutCol::RRAD]=radiative.get_evolution_real(
			OUTER_RADIUS
	).begin();
	__real_iterators[OutCol::MRAD]=radiative.get_evolution_real(
			OUTER_MASS
	).begin();
	__real_iterators[OutCol::ICONV_DERIV]=convective.get_evolution_real(
			MOMENT_OF_INERTIA_FIRST_DERIV
	).begin();
	__real_iterators[OutCol::IRAD_DERIV]=radiative.get_evolution_real(
			MOMENT_OF_INERTIA_FIRST_DERIV
	).begin();
	__real_iterators[OutCol::RSTAR_DERIV]=convective.get_evolution_real(
			OUTER_RADIUS_FIRST_DERIV
	).begin();
	__real_iterators[OutCol::RRAD_DERIV]=radiative.get_evolution_real(
			OUTER_RADIUS_FIRST_DERIV
	).begin();
	__real_iterators[OutCol::MRAD_DERIV]=radiative.get_evolution_real(
			OUTER_MASS_DERIV
	).begin();
	__real_iterators[OutCol::ICONV_SECOND_DERIV]=
		convective.get_evolution_real(
				MOMENT_OF_INERTIA_SECOND_DERIV
		).begin();
	__real_iterators[OutCol::IRAD_SECOND_DERIV]=
		radiative.get_evolution_real(
				MOMENT_OF_INERTIA_SECOND_DERIV
		).begin();
	__real_iterators[OutCol::RRAD_SECOND_DERIV]=radiative.get_evolution_real(
			OUTER_RADIUS_SECOND_DERIV
	).begin();
	__real_iterators[OutCol::SEMIMAJOR]=system.semimajor_evolution().begin();
	__real_iterators[OutCol::CONV_INCLINATION]=convective.get_evolution_real(
			INCLINATION
	).begin();
	__real_iterators[OutCol::RAD_INCLINATION]=radiative.get_evolution_real(
			INCLINATION
	).begin();
	__real_iterators[OutCol::CONV_PERIAPSIS]=convective.get_evolution_real(
			PERIAPSIS
	).begin();
	__real_iterators[OutCol::RAD_PERIAPSIS]=radiative.get_evolution_real(
			PERIAPSIS
	).begin();
	__real_iterators[OutCol::LCONV]=convective.get_evolution_real(
			ANGULAR_MOMENTUM
	).begin();
		__real_iterators[OutCol::LRAD]=radiative.get_evolution_real(
			ANGULAR_MOMENTUM
	).begin();
	if(__real_iterators[OutCol::AGE]==__last_age)
		fix_no_evolution(start_age, end_age, timestep, required_ages);
}

const ConstSolutionIterator &ConstSolutionIterator::operator++()
{
	for(unsigned i=0; i<__real_iterators.size(); ++i) 
		if(__real_iterators[i]!=__placeholder_iterator)
			++__real_iterators[i];
	++__mode;
	++__wind_saturation;
	ZoneOrientation 
		convective_orientation(*__real_iterators[OutCol::CONV_INCLINATION],
							   *__real_iterators[OutCol::CONV_PERIAPSIS]),
		radiative_orientation(*__real_iterators[OutCol::CONV_INCLINATION],
							  *__real_iterators[OutCol::CONV_PERIAPSIS]);
	__stellar_angmom=zone_to_zone_transform(
			radiative_orientation,
			convective_orientation,
			Eigen::Vector3d(0, 0, *__real_iterators[OutCol::LRAD])
	);
	__stellar_angmom[2]+=*__real_iterators[OutCol::LCONV];
	return *this;
}

double ConstSolutionIterator::real_quantity(OutCol::OutputColumns quantity)
{
	switch(quantity) {
		case OutCol::ITOT : 
			return (*__real_iterators[OutCol::ICONV]
					+
					*__real_iterators[OutCol::IRAD]);
		case OutCol::LSTAR : 
			return __star.luminosity(*__real_iterators[OutCol::AGE]);
		case OutCol::ITOT_DERIV :
			return (*__real_iterators[OutCol::ICONV_DERIV]
					+
					*__real_iterators[OutCol::IRAD_DERIV]);
		case OutCol::ITOT_SECOND_DERIV :
			return (*__real_iterators[OutCol::ICONV_SECOND_DERIV]
					+
					*__real_iterators[OutCol::IRAD_SECOND_DERIV]);
		case OutCol::WORB : return orbital_angular_velocity(
									__mstar,
									__mplanet,
									*__real_iterators[OutCol::SEMIMAJOR]
							);
		case OutCol::PORB : return 2.0*M_PI/real_quantity(OutCol::WORB);
		case OutCol::LORB : return orbital_angular_momentum(
									__mstar,
									__mplanet,
									*__real_iterators[OutCol::SEMIMAJOR],
									*__real_iterators[OutCol::ECCENTRICITY]
							);
		case OutCol::LTOT : return __stellar_angmom.norm();
		case OutCol::WSURF : return (*__real_iterators[OutCol::LCONV]
									 /
									 *__real_iterators[OutCol::ICONV]);
		case OutCol::WRAD : return (*__real_iterators[OutCol::LRAD]
									/
									*__real_iterators[OutCol::IRAD]);
		case OutCol::PSURF : return 2.0*M_PI/real_quantity(OutCol::WSURF);
		case OutCol::PRAD : return 2.0*M_PI/real_quantity(OutCol::WRAD);
		default : return *__real_iterators[quantity];
	}
}
