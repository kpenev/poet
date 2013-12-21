/**\file
 *
 * \brief Defines some of the methods of the classes for interpolating among
 * stellar evolution tracks.
 * 
 * \ingroup StellarSystem_group
 */

#include "StellarEvolution.h"
#include "Error.h"
#include <math.h>
#include <memory>
#include <sstream>

double LogArgDerivatives::transform_log_arg_deriv(unsigned order) const
{
	if(order==1) return underlying_deriv_values[0]/x;
	else if(order==2)
		return (underlying_deriv_values[1]-underlying_deriv_values[0])/(x*x);
	else throw Error::BadFunctionArguments(
			"Transforming log-derivatives of order higher than 2 is not "
			"implemented.");
}

double LogArgDerivatives::order(unsigned deriv_order) const
{
	if(deriv_order==0) {
		if(std::isnan(value)) value=calc_deriv(0);
		return value;
	}
	if(deriv_values.size()<deriv_order)
		for(unsigned i=deriv_values.size(); i<deriv_order; ++i)
			if(correct_log_arg) {
				underlying_deriv_values.push_back(calc_deriv(i+1));
				deriv_values.push_back(transform_log_arg_deriv(i+1));
			} else deriv_values.push_back(calc_deriv(i+1));
	return deriv_values[deriv_order-1];
}

double InterpolatedDerivatives::calc_deriv(unsigned deriv_order) const
{
	if(deriv_order>2) return 0.0;
	std::valarray<double> interp_values(interp_deriv->size());
	for(unsigned i=0; i<interp_deriv->size(); ++i)
		interp_values[i]=(*interp_deriv)[i]->order(deriv_order);
	return InterpolatingFunctionALGLIB(*interp_masses,
			interp_values)(stellar_mass);
}

double ScaledDerivatives::calc_deriv(unsigned deriv_order=1) const
{
	return (correct_log_arg ? 1.0 : 
			std::pow(scaling, static_cast<int>(deriv_order)))*
		underlying_deriv->order(deriv_order);
}

void EvolvingStellarQuantity::init_high_mass(
		const std::valarray<double> &masses_of_tracks,
		const std::list<const OneArgumentDiffFunction *> 
		&evolution_tracks)
{
	std::list<const OneArgumentDiffFunction *>::const_iterator track=
			evolution_tracks.begin();
	closest_high_mass=masses_of_tracks[0];
	double current_diff=std::abs(closest_high_mass-stellar_mass);
	for(size_t track_m_ind=0; track_m_ind<masses_of_tracks.size(); 
			++track_m_ind) {
		if(track==evolution_tracks.end())
			throw Error::BadFunctionArguments(
				" The number of evolution tracks provided "
				"is less than the number of masses in "
				"EvolvingStellarQuantity::init_high_mass.");
		double track_m=masses_of_tracks[track_m_ind],
		       diff=std::abs(track_m-stellar_mass);
		if(diff<=current_diff) {
			closest_high_mass=track_m;
			closest_high_mass_track=*track;
			current_diff=diff;
		}
		++track;
	}
	double age_scaling=std::pow(closest_high_mass/stellar_mass, age_scaling_high_mass);
	min_age=(initially_zero ? 0 : closest_high_mass_track->range_low());
	max_age=closest_high_mass_track->range_high();
	if(use_log_age) {
		if(!initially_zero) min_age=std::exp(min_age);
		max_age=std::exp(max_age);
	}
	min_age*=age_scaling;
	max_age*=age_scaling;
}

void EvolvingStellarQuantity::init_low_mass(
		const std::valarray<double> &masses_of_tracks,
		const std::list<const OneArgumentDiffFunction *> 
		&evolution_tracks,
		double max_low_mass)
{
	std::list<const OneArgumentDiffFunction *>::const_iterator 
		track=evolution_tracks.begin(),
		track_below=evolution_tracks.begin(),
		track_above=evolution_tracks.end();
	double mass_below=masses_of_tracks[0], mass_above=NaN;
	for(size_t i=0; i<masses_of_tracks.size(); ++i) {
		if(track==evolution_tracks.end())
			throw Error::BadFunctionArguments(
				"The number of evolution tracks provided "
				"is less than the number of masses in "
				"EvolvingStellarQuantity::init_low_mass.");
		double mass=masses_of_tracks[i];
		if(stellar_mass>=mass) {
			track_below=track;
			mass_below=mass;
		}
		if(stellar_mass<=mass && track_above==evolution_tracks.end()) {
			track_above=track;
			mass_above=mass;
		}
		low_mass_tracks.push_back(*(track++));
		low_masses.push_back(mass);
		if(mass>max_low_mass) break;
	}
	double age_scaling_below, age_scaling_above;
	if(masses_of_tracks.size()==1)
		age_scaling_below=age_scaling_above=0;
	else if(std::isnan(mass_above)) {
		std::ostringstream msg;
		msg << "Stellar mass: " << stellar_mass
			<< " is larger than the masses of all tracks.";
		throw Error::BadFunctionArguments(msg.str());
	} else {
		age_scaling_below=std::pow(mass_below/stellar_mass,
				age_scaling_low_mass);
		age_scaling_above=std::pow(mass_above/stellar_mass,
				age_scaling_low_mass);
	}
	if(track_above==evolution_tracks.end()) track_above--;
	double min_age_below=(*track_below)->range_low(),
		   min_age_above=(*track_above)->range_low(),
		   max_age_below=(*track_below)->range_high(),
		   max_age_above=(*track_above)->range_high();
	if(use_log_age) {
		min_age_below=std::exp(min_age_below);
		max_age_below=std::exp(max_age_below);
		min_age_above=std::exp(min_age_above);
		max_age_above=std::exp(max_age_above);
	}
	min_age=(initially_zero ? 0.0 : 
			std::max(age_scaling_below*min_age_below,
				age_scaling_above*min_age_above));
	max_age=std::min(age_scaling_below*max_age_below,
			age_scaling_above*max_age_above);
}

double EvolvingStellarQuantity::high_mass_interp(double age,
		const ScaledDerivatives **derivatives) const
{
	double age_scaling=std::pow(stellar_mass/closest_high_mass,
			age_scaling_high_mass), model_age=age*age_scaling;
	if(use_log_age) model_age=std::log(model_age);
	if(model_age<closest_high_mass_track->range_low()) {
		if(initially_zero) {
			if(derivatives!=NULL) *derivatives=NULL;
			return 0.0;
		} else {
			std::stringstream ss;
			ss << "Asking for out of range age " << age <<
					" in high_mass_interp";
			throw Error::BadFunctionArguments(ss.str());
			return NaN;
		}
	}
	if(model_age>closest_high_mass_track->range_high()*extrapolate_high_mass)
		throw Error::BadFunctionArguments("Asking for out of range age "
				"in high_mass_interp.");
	if(derivatives!=NULL) {
		*derivatives=new ScaledDerivatives(
				closest_high_mass_track->deriv(model_age), age_scaling,
				(use_log_age ? age : NaN) ,true);
		return (*derivatives)->order(0);
	} else return (*closest_high_mass_track)(model_age);
}

double EvolvingStellarQuantity::evaluate_track(double age,
		const OneArgumentDiffFunction &track,
		const FunctionDerivatives **derivatives) const
{
	bool too_young=age<track.range_low(),
		 too_old=age>track.range_high()*extrapolate_low_mass;
	if((!too_young || initially_zero) && !too_old) {
		if(derivatives==NULL) return (too_young ? 0.0 : track(age));
		else {
			if(too_young) *derivatives=new ZeroDerivatives;
			else if(use_log_age) *derivatives=new RemoveLogDeriv(
					std::exp(age), track.deriv(age), true);
			else *derivatives=track.deriv(age);
			return (*derivatives)->order(0);
		}
	} else return NaN;
}

double EvolvingStellarQuantity::low_mass_interp(double age,
		const FunctionDerivatives **derivatives) const
{
	int num_good_tracks=0;
	std::list<const OneArgumentDiffFunction *>::const_iterator 
			track=low_mass_tracks.begin();
	std::list<double>::const_iterator mass=low_masses.begin();
	std::list<double> scaled_ages;
	for(; track!=low_mass_tracks.end(); ++track) {
		if(*mass==stellar_mass || low_mass_tracks.size()==1)
			return evaluate_track((use_log_age ? std::log(age) : age),
					**track, derivatives);
		double model_age=std::pow(stellar_mass/(*mass), 
				age_scaling_low_mass)*age;
		if(use_log_age) model_age=std::log(model_age);
		scaled_ages.push_back(model_age);
		if(((*track)->range_low()<=model_age || initially_zero) && 
				(*track)->range_high()*extrapolate_low_mass>=model_age)
			++num_good_tracks;
		++mass;
	}
	std::valarray<double> 
		*good_masses=new std::valarray<double>(num_good_tracks), 
		good_values(num_good_tracks);
	std::valarray<const FunctionDerivatives *> *good_derivatives=
		new std::valarray<const FunctionDerivatives *>(num_good_tracks);
	mass=low_masses.begin();

	std::list<double>::const_iterator model_age=scaled_ages.begin();
	int good_ind=0;
	for(track=low_mass_tracks.begin(); track!=low_mass_tracks.end(); 
			++track) {
		double value=evaluate_track(*model_age, **track,
				(derivatives ? &((*good_derivatives)[good_ind]) : NULL));
		if(!std::isnan(value)) {
			good_values[good_ind]=value;
			(*good_masses)[good_ind]=*mass;
			++good_ind;
		}
		++mass;
		++model_age;
	}
	double result;
	if(derivatives==NULL) {
		result=InterpolatingFunctionALGLIB(*good_masses, good_values)(
				stellar_mass);
		delete good_masses;
		delete good_derivatives;
	} else {
		*derivatives=new InterpolatedDerivatives(stellar_mass,
				good_derivatives, good_masses, (use_log_age ? age : NaN),
				true);
		result=(*derivatives)->order(0);
	}
	return result;
}

EvolvingStellarQuantity::EvolvingStellarQuantity(
		double mass, const std::valarray<double> &masses_of_tracks,
		const std::list<const OneArgumentDiffFunction *> &evolution_tracks,
		bool log_age, double max_low_mass, double low_mass_age_scaling,
		double high_mass_age_scaling, double low_mass_extrapolate,
		double high_mass_extrapolate, bool starts_zero) :
	use_log_age(log_age), initially_zero(starts_zero), stellar_mass(mass), 
	age_scaling_low_mass(low_mass_age_scaling),
	age_scaling_high_mass(high_mass_age_scaling),
	extrapolate_low_mass(low_mass_extrapolate),
	extrapolate_high_mass(high_mass_extrapolate), 
	is_low_mass(mass<=max_low_mass)
{
	using namespace std;
	if(is_low_mass) {
		init_low_mass(masses_of_tracks, evolution_tracks, max_low_mass);
	}
	else {
		init_high_mass(masses_of_tracks, evolution_tracks);
	}
}

double EvolvingStellarQuantity::operator()(double age) const
{
	if(is_low_mass) return low_mass_interp(age);
	else return high_mass_interp(age);
}

const FunctionDerivatives *EvolvingStellarQuantity::deriv(double age) const
{
	if(is_low_mass) {
		const FunctionDerivatives *deriv;
		low_mass_interp(age, &deriv);
		return deriv;
	} else {
		const ScaledDerivatives *deriv;
		high_mass_interp(age, &deriv);
		if(deriv==NULL) return new ZeroDerivatives;
		return deriv;
	}
}

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
				*Iconv_iter, std::valarray<double>(), smooth_conv_inertia,
				conv_inertia_nodes));
		std::cout << "done" << std::endl;

		size_t first_core_index=0;
		while((*Mrad_iter)[first_core_index+1]==0) ++first_core_index;
		std::slice core_slice(first_core_index,
				log_ages.size()-first_core_index, 1);
		std::valarray<double> core_log_ages=log_ages[core_slice];
		core_formation=(*ages_iter)[first_core_index];
		std::cout << "core_formation=" << core_formation << std::endl;

		std::cout << "Interpolating Irad...";
		std::cout.flush();
		interpolated_rad_inertia.push_back(
			new InterpolatingFunctionALGLIB(core_log_ages,
				(*Irad_iter)[core_slice], std::valarray<double>(),
				smooth_rad_inertia/tabulated_masses[i], rad_inertia_nodes));
		std::cout << "done" << std::endl;

		std::cout << "Interpolating Mrad...";
		std::cout.flush();
		interpolated_rad_mass.push_back(
			new InterpolatingFunctionALGLIB(core_log_ages,
				(*Mrad_iter)[core_slice], std::valarray<double>(),
				smooth_rad_mass/tabulated_masses[i], rad_mass_nodes));
		std::cout << "done" << std::endl;
		std::cout << "Interpolating Rrad...";
		std::cout.flush();
		interpolated_core_env_boundary.push_back(
			new InterpolatingFunctionALGLIB(core_log_ages, 
				(*core_boundary_iter)[core_slice], std::valarray<double>(),
				smooth_core_env_boundary, core_env_boundary_nodes));
		std::cout << "done" << std::endl;
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
		"evolution quantity.");
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
