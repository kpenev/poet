/**\file
 *
 * \brief Defines some of the methods of the classes for generating stellar
 * evolution interpolators from the YREC tracks. 
 * 
 * \ingroup StellarSystem_group
 */

#include "YRECIO.h"
#include <string>
#include <sstream>
#include <algorithm>

const std::streamsize max_line_length=1000;

EvolutionIterator::EvolutionIterator(const EvolutionIterator &orig) :
	mass_iter(orig.mass_iter), age_iter(orig.age_iter),
	radius_iter(orig.radius_iter), luminosity_iter(orig.luminosity_iter),
	rad_mass_iter(orig.rad_mass_iter),
	core_boundary_iter(orig.core_boundary_iter),
	conv_inertia_iter(orig.conv_inertia_iter),
	rad_inertia_iter(orig.rad_inertia_iter) {}

EvolutionIterator &EvolutionIterator::operator=(const EvolutionIterator &rhs)
{
	mass_iter=rhs.mass_iter;
	age_iter=rhs.age_iter;
	radius_iter=rhs.radius_iter;
	luminosity_iter=rhs.luminosity_iter;
   	rad_mass_iter=rhs.rad_mass_iter; 
	core_boundary_iter=rhs.core_boundary_iter;
	conv_inertia_iter=rhs.conv_inertia_iter;
	rad_inertia_iter=rhs.rad_inertia_iter;
	return *this;
}

EvolutionIterator &EvolutionIterator::operator++()
{
	++mass_iter; ++age_iter; ++radius_iter; ++luminosity_iter;
	++rad_mass_iter; ++core_boundary_iter; ++conv_inertia_iter;
	++rad_inertia_iter;
	return *this;
}

EvolutionIterator EvolutionIterator::operator++(int)
{
	EvolutionIterator result(*this); ++(*this); return result;
}

bool EvolutionIterator::operator==(const EvolutionIterator &rhs)
{
	return mass_iter==rhs.mass_iter;
}

YRECHeader::YRECHeader(std::ifstream &track, const std::string &filename)
{
	const std::string mass_word="Mtot/Msun", age_word="Age(Gyr)", 
	      radius_word="log(R/Rsun)", log_luminosity_word="log(L/Lsun)",
		  env_mass_word="m_envp/M", env_radius_word1="r_envp/M",
		  env_radius_word2="r_envp/R", Irad_word="I_rad", Iconv_word="I_env";
	char word[max_line_length+1];
	char line_characters[max_line_length+1];
	while(track.get()=='#') {
		track.getline(line_characters, max_line_length);
		if(!track.good()) 
			throw Error::IO(filename, "while reading header.");
		std::istringstream line(line_characters);
		assert(line.good());
		for(int word_index=0; !line.eof(); word_index++) {
			assert(line.good());
			line >> word;
			if(word==mass_word) {
				std::string error_msg="failed to read track mass: ";
				if(!line.good())
					throw Error::IO(filename, error_msg +
							(line.eof() ? "end of line before =" :
							 "read error before ="));
				line >> word;
				if(!line.good()) throw Error::IO(filename, error_msg +
						(line.eof() ? "end of line after =" :
						 "read error at ="));
				line >> track_mass;
				if(!line.good()) throw Error::IO(filename, error_msg +
						(line.eof() ? "end of line before value" :
						"error parsing value"));
			}
			if(word==age_word) age_col=word_index;
			else if(word==radius_word) radius_col=word_index;
			else if(word==log_luminosity_word) log_luminosity_col=word_index;
			else if(word==env_mass_word) 
				envelope_mass_col=word_index;
			else if(word==env_radius_word1 || word==env_radius_word2)
				rad_conv_boundary_col=word_index;
			else if(word==Irad_word) rad_inertia_col=word_index;
			else if(word==Iconv_word) conv_inertia_col=word_index;
		}
	}
	track.unget();
}

void YRECEvolution::read_model_file(const std::string &filename)
{
	std::cerr << "Reading: " << filename << std::endl;
	std::ifstream file(filename.c_str());
	YRECHeader header(file, filename);
	mass_list.push_back(header.get_mass());
	int required_columns=std::max(header.get_age_col(), 
			header.get_radius_col());
	required_columns=std::max(required_columns,
			header.get_log_luminosity_col());
	required_columns=std::max(required_columns,
			header.get_envelope_mass_col());
	required_columns=std::max(required_columns, 
			header.get_core_boundary_col());
	required_columns=std::max(required_columns,
			header.get_rad_inertia_col());
	required_columns=std::max(required_columns,
			header.get_conv_inertia_col());
	char line_characters[max_line_length+1];
	std::list<double> track_ages, track_radii, track_luminosities,
		track_rad_masses, track_core_boundaries, track_conv_inertias,
		track_rad_inertias;
	double inertia_norm=AstroConst::solar_mass*std::pow(
			AstroConst::solar_radius, 2)*1e7;
	for(int line_number=0; !file.eof(); line_number++) {
		std::cerr << "Reading line " << line_number << "\r";
		std::cerr.flush();
		std::ostringstream error_message;
		error_message << "while reading data line " << line_number;
		if(!file.good()) throw Error::IO(filename, error_message.str());
		file.getline(line_characters, max_line_length);
		if(file.gcount()==0) continue;
		std::istringstream line(line_characters);
		int column=0;
		for(; !line.eof(); column++) {
			double value;
			line >> value;
			if(column==header.get_age_col()) track_ages.push_back(value);
			else if(column==header.get_radius_col())
				track_radii.push_back(std::pow(10.0, value));
			else if(column==header.get_log_luminosity_col())
				track_luminosities.push_back(std::pow(10.0, value));
			else if(column==header.get_envelope_mass_col())
				track_rad_masses.push_back(value*header.get_mass());
			else if(column==header.get_core_boundary_col())
				track_core_boundaries.push_back(value*track_radii.back());
			else if(column==header.get_rad_inertia_col()) 
				track_rad_inertias.push_back(value/inertia_norm);
			else if(column==header.get_conv_inertia_col())
				track_conv_inertias.push_back(value/inertia_norm);
		}
		if(column<=required_columns) 
			throw Error::IO(filename, error_message.str());
	}
	ages.push_back(list_to_valarray(track_ages));
	radii.push_back(list_to_valarray(track_radii));
	luminosities.push_back(list_to_valarray(track_luminosities));
	rad_masses.push_back(list_to_valarray(track_rad_masses));
	core_boundaries.push_back(list_to_valarray(track_core_boundaries));
	rad_inertias.push_back(list_to_valarray(track_rad_inertias));
	conv_inertias.push_back(list_to_valarray(track_conv_inertias));
}

EvolutionIterator YRECEvolution::begin()
{
	EvolutionIterator result;
	result.mass_iter=mass_list.begin();
	result.age_iter=ages.begin();
	result.radius_iter=radii.begin();
	result.luminosity_iter=luminosities.begin();
	result.rad_mass_iter=rad_masses.begin();
	result.core_boundary_iter=core_boundaries.begin();
	result.conv_inertia_iter=conv_inertias.begin();
	result.rad_inertia_iter=rad_inertias.begin();
	return result;
}

EvolutionIterator YRECEvolution::end()
{
	EvolutionIterator result;
	result.mass_iter=mass_list.end();
	result.age_iter=ages.end();
	result.radius_iter=radii.end();
	result.luminosity_iter=luminosities.end();
	result.rad_mass_iter=rad_masses.end();
	result.core_boundary_iter=core_boundaries.end();
	result.conv_inertia_iter=conv_inertias.end();
	result.rad_inertia_iter=rad_inertias.end();
	return result;
}

void YRECEvolution::move(EvolutionIterator &dest, EvolutionIterator &source)
{
	mass_list.splice(dest.mass_iter, mass_list, source.mass_iter);
	ages.splice(dest.age_iter, ages, source.age_iter);
	radii.splice(dest.radius_iter, radii, source.radius_iter);
	luminosities.splice(dest.luminosity_iter, luminosities,
			source.luminosity_iter);
	rad_masses.splice(dest.rad_mass_iter, rad_masses, source.rad_mass_iter);
	core_boundaries.splice(dest.core_boundary_iter, core_boundaries,
			source.core_boundary_iter);
	conv_inertias.splice(dest.conv_inertia_iter, conv_inertias,
			source.conv_inertia_iter);
	rad_inertias.splice(dest.rad_inertia_iter, rad_inertias,
			source.rad_inertia_iter);
}

void YRECEvolution::sort_masses()
{
	EvolutionIterator iter=begin(), stop_iter=end();
	double largest_sorted=*iter.mass_iter;
	while(iter!=stop_iter) {
		while(*iter.mass_iter>=largest_sorted && iter!=stop_iter) {
			largest_sorted=*iter.mass_iter;
			++iter;
		}
		if(iter==stop_iter) break;
		EvolutionIterator dest=begin();
		while(*(dest.mass_iter)<=*iter.mass_iter) {++dest;}
		EvolutionIterator source=iter++;
		move(dest, source);
	}
}

YRECEvolution::YRECEvolution(const std::string &model_directory,
	double smooth_conv_inertia, double smooth_rad_inertia,
	double smooth_rad_mass, int conv_inertia_nodes,
	int rad_inertia_nodes, int rad_mass_nodes)
{
	DIR *dirstream=opendir(model_directory.c_str());
	std::string join;
	if(model_directory[model_directory.size()-1]=='/') join="";
	else join="/";
	if(dirstream==NULL) throw Error::PathNotFound(
		"in YRECEvolution constructor.", model_directory);
	struct dirent *entry; 
	while((entry=readdir(dirstream))) {
		std::string fname(entry->d_name);
		if(fname[0]!='.' && fname.substr(fname.size()-6)==".track")
			read_model_file(model_directory+join+entry->d_name);
	}
	if(closedir(dirstream)) throw Error::Runtime(
		"Failed to close directory stream tied to "+model_directory+
		" in YRECEvolution constructor.");
	sort_masses();
	std::valarray<double> masses = list_to_valarray(mass_list);
	interpolate_from(masses, ages, radii, conv_inertias, rad_inertias,
			rad_masses, core_boundaries, NaN, smooth_conv_inertia,
			smooth_rad_inertia, smooth_rad_mass, NaN, 0, conv_inertia_nodes,
			rad_inertia_nodes, rad_mass_nodes, 0, luminosities, NaN, 0);
}
