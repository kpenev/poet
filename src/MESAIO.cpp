#include "MESAIO.h"

std::ostream &operator<<(std::ostream &os, MESA::Column col)
{
	switch(col) {
		case MESA::AGE : os << "AGE"; break;
		case MESA::LOG_RSTAR : os << "LOG_RSTAR"; break;
		case MESA::LOG_LSTAR : os << "LOG_LSTAR"; break;
		case MESA::MRAD : os << "MRAD"; break;
		case MESA::RRAD : os << "RRAD"; break;
		case MESA::ICONV : os << "ICONV"; break;
		case MESA::IRAD : os << "IRAD"; break;
		default : throw Error::IO("Unrecognized MESA column encountered!");
	};
	return os;
}

namespace MESA {

	void Header::set_column_names()
	{
		__column_names.resize(NUM_COLUMNS);
		__column_names[AGE]="star_age";
		__column_names[LOG_RSTAR]="log_R";
		__column_names[LOG_LSTAR]="log_L";
		__column_names[MRAD]="conv_mx1_bot";
		__column_names[RRAD]="conv_mx1_bot_r";
		__column_names[ICONV]="Irot_conv";
		__column_names[IRAD]="Irot_rad";
	}

	void Header::read_column_numbers(std::istream &track,
			const std::string &filename, unsigned &line_number)
	{
		std::string line="";
		for(;line==""; ++line_number) {
			if(!track.good()) {
				std::ostringstream msg;
				msg << "Failed to extract line " << line_number
					<< " from " << filename << ".";
				throw Error::IO(msg.str());
			}
			std::getline(track, line);
		}
		std::istringstream line_parser(line);
		int input_int;
		for(int word_ind=1; line_parser.good(); ++word_ind) {
			line_parser >> input_int;
			if(line_parser.good() && input_int!=word_ind) {
				std::ostringstream message;
				message << "The contents of line " << line_number
					<< " of " << filename << " are not sequential integers "
					"starting with 1.";
				throw Error::IO(message.str());
			}
		}
	}

	Header::Header(std::ifstream &track, const std::string &filename)
	{
		set_column_names();
		unsigned line_number=0;
		read_column_numbers(track, filename, line_number);
		std::string line;

		++line_number;
		std::getline(track, line);
		std::istringstream line_parser(line);
		int mass_column;
		const std::string mass_column_name="initial_mass";
		std::string word;
		for(mass_column=0; line_parser.good() && word!=mass_column_name;
				++mass_column)
			line_parser >> word;
		if(!line_parser.good()) {
			std::ostringstream message;
			message << "Line " << line_number << " of " << filename 
				<< " does not contain the word '" << mass_column_name << "'";
			throw Error::IO(message.str());
		}

		++line_number;
		std::getline(track, line);

		for(line_parser.str(line); line_parser.good() && mass_column>0;
				--mass_column) line_parser >> __track_mass;
		if(!line_parser.good()) {
			std::ostringstream message;
			message << "Failed parsing the model mass from line "
				<< line_number << " of " + filename + ".";
			throw Error::IO(message.str());
		}

		read_column_numbers(track, filename, line_number);

		++line_number;
		std::getline(track, line);
		line_parser.str(line);
		__column_numbers.resize(NUM_COLUMNS, -1);
		for(int column=0; true; ++column) {
			line_parser >> word;
			if(!line_parser.good()) break;
			for(int i=0; i<NUM_COLUMNS; ++i)
				if(word==__column_names[i])
					__column_numbers[i]=column;
		}
		for(int i=0; i<NUM_COLUMNS; ++i)
			if(__column_numbers[i]==-1) {
				std::ostringstream message;
				message << "Failed to find '" << __column_names[i]
					<< "' on line " << line_number << " of " << filename
					<< ".";
				throw Error::IO(message.str());
			}
	}
	EvolutionIterator &EvolutionIterator::operator=(
			const EvolutionIterator &rhs)
	{
		mass_iter=rhs.mass_iter;
		quantity_iter=rhs.quantity_iter;
		return *this;
	}

	EvolutionIterator &EvolutionIterator::operator++()
	{
		++mass_iter;
		for(size_t i=0; i<quantity_iter.size(); ++i) ++quantity_iter[i];
		return *this;
	}

	///Reads a single evolution track file
	void Evolution::read_model_file(const std::string &filename)
	{
		std::ifstream track(filename.c_str());
		Header header(track, filename);
		std::valarray< std::list<double> > track_columns=parse_columns(track,
				header.get_all_columns());
		__mass_list.push_back(header.get_mass());
		const double Inorm=AstroConst::solar_mass*
			std::pow(AstroConst::solar_radius, 2);
		for(int i=0; i<NUM_COLUMNS; ++i)
			if(i==AGE)
				track_quantities[i].push_back(
						list_to_valarray(track_columns[i])*1e-9);
			else if(i==LOG_RSTAR || i==LOG_LSTAR)
				track_quantities[i].push_back(
						std::pow(10.0, list_to_valarray(track_columns[i])));
			else if(i==ICONV || i==IRAD)
				track_quantities[i].push_back(
						list_to_valarray(track_columns[i])/Inorm);
			else track_quantities[i].push_back(
						list_to_valarray(track_columns[i]));
	}

	EvolutionIterator Evolution::begin()
	{
		EvolutionIterator result;
		result.mass_iter=__mass_list.begin();
		for(size_t i=0; i<track_quantities.size(); ++i)
			result.quantity_iter[i]=track_quantities[i].begin();
		return result;
	}

	EvolutionIterator Evolution::end()
	{
		EvolutionIterator result;
		result.mass_iter=__mass_list.end();
		for(size_t i=0; i<track_quantities.size(); ++i)
			result.quantity_iter[i]=track_quantities[i].end();
		return result;
	}

	void Evolution::move(EvolutionIterator &dest, EvolutionIterator &source)
	{
		__mass_list.splice(dest.mass_iter, __mass_list, source.mass_iter);
		for(size_t i=0; i<track_quantities.size(); ++i)
			track_quantities[i].splice(dest.quantity_iter[i],
					track_quantities[i], source.quantity_iter[i]);
	}

	///Sorts the data by mass.
	void Evolution::sort_masses()
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

	Evolution::Evolution(const std::string &model_directory,
			double smooth_conv_inertia, double smooth_rad_inertia,
			double smooth_rad_mass, int conv_inertia_nodes,
			int rad_inertia_nodes, int rad_mass_nodes) :
		track_quantities(NUM_COLUMNS)
	{
		DIR *dirstream=opendir(model_directory.c_str());
		std::string join;
		if(model_directory[model_directory.size()-1]=='/') join="";
		else join="/";
		if(dirstream==NULL) throw Error::PathNotFound(
				"in MESA::Evolution constructor.", model_directory);
		struct dirent *entry; 
		while((entry=readdir(dirstream))) {
			std::string fname(entry->d_name);
			if(fname[0]!='.' && fname.substr(fname.size()-5)==".data")
				read_model_file(model_directory+join+entry->d_name);
		}
		if(closedir(dirstream)) throw Error::Runtime(
				"Failed to close directory stream tied to "+model_directory+
				" in MESA::Evolution constructor.");
		sort_masses();
		std::valarray<double> masses = list_to_valarray(__mass_list);
		interpolate_from(masses, track_quantities[AGE],
				track_quantities[RSTAR], track_quantities[ICONV],
				track_quantities[IRAD], track_quantities[MRAD],
				track_quantities[RRAD], NaN, smooth_conv_inertia,
				smooth_rad_inertia, smooth_rad_mass, NaN, 0,
				conv_inertia_nodes, rad_inertia_nodes, rad_mass_nodes, 0,
				track_quantities[LSTAR], NaN, 0);
	}
};
