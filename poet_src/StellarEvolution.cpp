/**\file
 *
 * \brief Defines some of the methods of the StellarEvolution class.
 * 
 * \ingroup StellarSystem_group
 */

void StellarEvolution::interpolate_from(
	const std::valarray<double> &tabulated_masses,
	const std::list< std::valarray<double> > &tabulated_ages,
	const std::list< std::valarray<double> > &tabulated_radii,
	const std::list< std::valarray<double> > &tabulated_conv_inertia,
	const std::list< std::valarray<double> > &tabulated_rad_inertia,
	const std::list< std::valarray<double> > &tabulated_rad_mass,
	const std::list< std::valarray<double> > &tabulated_core_env_boundary,

	double smooth_radius,
	double smooth_conv_inertia,
	double smooth_rad_inertia,
	double smooth_rad_mass,
	double smooth_core_env_boundary,

	int radius_nodes,
	int conv_inertia_nodes,
	int rad_inertia_nodes,
	int rad_mass_nodes,
	int core_env_boundary_nodes,

	const std::list< std::valarray<double> > &tabulated_luminosities,
	double smooth_luminosities,
	int luminosities_nodes,

	double max_low_mass, double low_mass_age_scaling,
	double high_mass_age_scaling, double low_mass_extrapolate,
	double high_mass_extrapolate)
{
	mass_break=max_low_mass;
	low_age_scaling=low_mass_age_scaling;
	high_age_scaling=high_mass_age_scaling;
	extrapolate_low=low_mass_extrapolate;
	extrapolate_high=high_mass_extrapolate;

	track_masses=new std::valarray<double>(tabulated_masses);
	size_t num_tracks=tabulated_masses.size();
	std::list< std::valarray<double> >::const_iterator 
		ages_iter=tabulated_ages.begin(),
		radii_iter=tabulated_radii.begin(),
		Iconv_iter=tabulated_conv_inertia.begin(),
		Irad_iter=tabulated_rad_inertia.begin(),
		Mrad_iter=tabulated_rad_mass.begin(),
		core_boundary_iter=tabulated_core_env_boundary.begin(),
		L_iter=tabulated_luminosities.begin();

	for(size_t i=0; i<num_tracks; ++i) {
		std::cout << "Track " << i << std::endl;
		if(ages_iter==tabulated_ages.end())
			throw Error::BadFunctionArguments(
				"The number of age arrays is smaller than the number of "
				"masses in StellarEvolution::interpolate_from.");
		if(radii_iter==tabulated_radii.end())
			throw Error::BadFunctionArguments(
				"The number of radii arrays is smaller than the number of "
				"masses in StellarEvolution::interpolate_from.");
		if(Iconv_iter==tabulated_conv_inertia.end())
			throw Error::BadFunctionArguments(
				"The number of convective moment of inertia arrays is "
				"smaller than the number of masses in "
				"StellarEvolution::interpolate_from.");
		if(Irad_iter==tabulated_rad_inertia.end())
			throw Error::BadFunctionArguments(
				"The number of radiative moment of inertia arrays is smaller"
				" than the number of masses in "
				"StellarEvolution::interpolate_from.");
		if(Mrad_iter==tabulated_rad_mass.end())
			throw Error::BadFunctionArguments(
				"The number of radiative mass arrays is smaller than the "
				"number of masses in StellarEvolution::interpolate_from.");
		if(core_boundary_iter==tabulated_core_env_boundary.end())
			throw Error::BadFunctionArguments(
				"The number of core-envelope boundary arrays is smaller than"
				" the number of masses in "
				"StellarEvolution::interpolate_from.");
		if(tabulated_luminosities.size()>0 &&
				L_iter==tabulated_luminosities.end())
			throw Error::BadFunctionArguments(
				"The number of luminosity arrays is smaller than the number "
				"of masses in StellarEvolution::interpolate_from.");
		std::valarray<double> log_ages=std::log(*ages_iter);

		std::cout << "Interpolating radius (num_ages: " << log_ages.size()
			<< ", num_radii: " << radii_iter->size() << ", smoothing: "
			<< smooth_radius << ", nodes: " << radius_nodes <<")...";
		std::cout.flush();
		interpolated_radius.push_back(
			new InterpolatingFunctionALGLIB(log_ages, *radii_iter,
				std::valarray<double>(), smooth_radius, radius_nodes));
		std::cout << "done" << std::endl;

		std::cout << "Interpolating Iconv...";
		std::cout.flush();
		interpolated_conv_inertia.push_back(
			new InterpolatingFunctionALGLIB(log_ages,
                                            *Iconv_iter,
                                            std::valarray<double>(),
                                            smooth_conv_inertia,
                                            conv_inertia_nodes)
        );
		std::cout << "done" << std::endl;

		size_t first_core_index=0;
		while(
				!tabulated_rad_mass.empty()
				&&
				first_core_index + 1 < Mrad_iter->size()
				&&
				(*Mrad_iter)[first_core_index+1]==0
		) ++first_core_index;
		if(
            Mrad_iter->size() == 0
            || 
            first_core_index + 1 == Mrad_iter->size() 
            ||
            (*Mrad_iter)[first_core_index + 1]==0
        ) {
			core_formation=Inf;
			interpolated_rad_inertia.push_back(new ZeroFunction());
			interpolated_rad_mass.push_back(new ZeroFunction());
			interpolated_core_env_boundary.push_back(new ZeroFunction());
		} else {
			std::slice core_slice(first_core_index,
					log_ages.size()-first_core_index, 1);
			std::valarray<double> core_log_ages=log_ages[core_slice];
			core_formation=(*ages_iter)[first_core_index];
			std::cout << "core_formation=" << core_formation << std::endl;

			std::cout << "Interpolating Irad...";
			std::cout.flush();
			interpolated_rad_inertia.push_back(
					new InterpolatingFunctionALGLIB(
						core_log_ages,
						(*Irad_iter)[core_slice],
						std::valarray<double>(),
						smooth_rad_inertia/tabulated_masses[i],
						rad_inertia_nodes
					)
			);
			std::cout << "done" << std::endl;

			std::cout << "Interpolating Mrad...";
			std::cout.flush();
			interpolated_rad_mass.push_back(
					new InterpolatingFunctionALGLIB(
						core_log_ages,
						(*Mrad_iter)[core_slice],
						std::valarray<double>(),
						smooth_rad_mass/tabulated_masses[i],
						rad_mass_nodes
					)
			);
			std::cout << "done" << std::endl;
			std::cout << "Interpolating Rrad...";
			std::cout.flush();
			interpolated_core_env_boundary.push_back(
					new InterpolatingFunctionALGLIB(
						core_log_ages, 
						(*core_boundary_iter)[core_slice],
						std::valarray<double>(),
						smooth_core_env_boundary,
						core_env_boundary_nodes
					)
			);
			std::cout << "done" << std::endl;
		}
		if(tabulated_luminosities.size()>0) {
			std::cout << "Interpolating Luminosity...";
			std::cout.flush();
			interpolated_luminosity.push_back(
					new InterpolatingFunctionALGLIB(log_ages, *L_iter,
						std::valarray<double>(), smooth_luminosities,
						luminosities_nodes));
			std::cout << "done" << std::endl;
			++L_iter;
		}
		++ages_iter;
		++radii_iter;
		++Iconv_iter;
		++Irad_iter;
		++Mrad_iter;
		++core_boundary_iter;

	}
	if(ages_iter!=tabulated_ages.end())
		throw Error::BadFunctionArguments(
			"The number of age arrays is larger than the number of masses in"
			" StellarEvolution::interpolate_from.");
	if(radii_iter!=tabulated_radii.end())
		throw Error::BadFunctionArguments(
			"The number of radii arrays is larger than the number of masses "
			"in StellarEvolution::interpolate_from.");
	if(Iconv_iter!=tabulated_conv_inertia.end())
		throw Error::BadFunctionArguments(
			"The number of convective moment of inertia arrays is larger "
			"than the number of masses in "
			"StellarEvolution::interpolate_from.");
	if(Irad_iter!=tabulated_rad_inertia.end())
		throw Error::BadFunctionArguments(
			"The number of radiative moment of inertia arrays is larger than"
			" the number of masses in StellarEvolution::interpolate_from.");
	if(Mrad_iter!=tabulated_rad_mass.end())
		throw Error::BadFunctionArguments(
			"The number of radiative mass arrays is larger than the number "
			"of masses in StellarEvolution::interpolate_from.");
	if(core_boundary_iter!=tabulated_core_env_boundary.end())
		throw Error::BadFunctionArguments(
			"The number of core-envelope boundary arrays is larger than the "
			"number of masses in StellarEvolution::interpolate_from.");
	if(L_iter!=tabulated_luminosities.end())
		throw Error::BadFunctionArguments(
			"The number of log luminosity arrays is larger than the number "
			"of masses in StellarEvolution::interpolate_from.");
}

const EvolvingStellarQuantity 
	*StellarEvolution::interpolate_moment_of_inertia(
			double stellar_mass, StellarZone zone, double) const
{
	switch (zone) {
		case radiative : return new EvolvingStellarQuantity(stellar_mass,
					     *track_masses, 
					     interpolated_rad_inertia, true, mass_break,
						 low_age_scaling, high_age_scaling, extrapolate_low,
						 extrapolate_high, true);
		case convective : return new EvolvingStellarQuantity(
						  stellar_mass,
						  *track_masses,
						  interpolated_conv_inertia, true, mass_break,
						  low_age_scaling, high_age_scaling, extrapolate_low,
						  extrapolate_high);
		case total :  return new SumQuantity(
								  new EvolvingStellarQuantity(stellar_mass,
									  *track_masses, 
									  interpolated_rad_inertia, true,
									  mass_break, low_age_scaling,
									  high_age_scaling, extrapolate_low,
									  extrapolate_high, true),
								  new EvolvingStellarQuantity(
									  stellar_mass,
									  *track_masses,
									  interpolated_conv_inertia, true,
									  mass_break, low_age_scaling,
									  high_age_scaling, extrapolate_low,
									  extrapolate_high)
								  );
		default: throw Error::BadFunctionArguments(
			 "Moment of inertia requested for an unrecognized "
			 "stellar zone.");
	}
}

const EvolvingStellarQuantity *StellarEvolution::interpolate_radius(
		double stellar_mass, double) const
{
	return new EvolvingStellarQuantity(stellar_mass, *track_masses, 
			interpolated_radius, true, mass_break, low_age_scaling,
			high_age_scaling, extrapolate_low, extrapolate_high);
}

const EvolvingStellarQuantity *StellarEvolution::interpolate_luminosity(
		double stellar_mass, double) const
{
	if(interpolated_luminosity.size()==0) return NULL;
	return new EvolvingStellarQuantity(stellar_mass, *track_masses, 
			interpolated_luminosity, true, mass_break, low_age_scaling,
			high_age_scaling, extrapolate_low, extrapolate_high);
}

const EvolvingStellarQuantity *StellarEvolution::interpolate_zone_mass(
		double stellar_mass, StellarZone zone) const
{
	if (zone!=radiative) throw Error::BadFunctionArguments(
		"Only the radiative zone mass is supported as a stellar "
		"evolution quantity."
    );
	return new EvolvingStellarQuantity(stellar_mass, *track_masses,
			interpolated_rad_mass, true, mass_break, low_age_scaling,
			high_age_scaling, extrapolate_low, extrapolate_high, true);
}

const EvolvingStellarQuantity *StellarEvolution::interpolate_core_boundary(
		double stellar_mass, double) const
{
	return new EvolvingStellarQuantity(stellar_mass, *track_masses,
			interpolated_core_env_boundary, true, mass_break,
			low_age_scaling, high_age_scaling, extrapolate_low,
			extrapolate_high, true);
}

#ifndef NO_SERIALIZE
void StellarEvolution::load_state(const std::string &filename)
{
    std::ifstream ifs(filename.c_str());
    boost::archive::text_iarchive ia(ifs);
    ia >> (*this);
    ifs.close();
}

void StellarEvolution::save_state(const std::string &filename) const
{
	std::ofstream ofs(filename.c_str());
	boost::archive::text_oarchive oa(ofs);
	oa << (*this);
}
#endif
