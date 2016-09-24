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
        __column_names[MTRACK]="M_ini";
        __column_names[AGE]="age";
        __column_names[LOG_RSTAR]="";
        __column_names[RSTAR]="R_star";
        __column_names[LOG_LSTAR]="";
        __column_names[LSTAR]="L_star";
        __column_names[MRAD]="M_rad";
        __column_names[RRAD]="R_tachocline";
        __column_names[ICONV]="I_conv";
        __column_names[IRAD]="I_rad";
    }

    void Header::read_column_numbers(std::istream &track,
                                     const std::string &filename,
                                     unsigned &line_number)
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

        std::cout << "Line: " << line << std::endl;
        std::istringstream line_parser(line);
        std::string column_name;
        for(
            int column_number = 0;
            line_parser.good();
            ++column_number
        ) {
            std::getline(line_parser, column_name, ',');
            int colname_index = 0;
            for(
                    ;
                    (
                        colname_index < NUM_COLUMNS
                        &&
                        column_name != __column_names[colname_index] 
                    );
                    ++colname_index
            ) {}
            if(colname_index < NUM_COLUMNS)
                __column_numbers[colname_index] = column_number;
        }

    }

    Header::Header(std::ifstream &track, const std::string &filename)
    {
        set_column_names();
        __column_numbers.resize(NUM_COLUMNS, -1);

        unsigned line_number=0;
        read_column_numbers(track, filename, line_number);

        for(int i=0; i<NUM_COLUMNS; ++i)
            if(
                __column_numbers[i]==-1
                &&
                !(
                    (i == LOG_RSTAR && __column_numbers[RSTAR] != -1) 
                    ||
                    (i == RSTAR && __column_numbers[LOG_RSTAR] != -1)
                    ||
                    (i == LOG_LSTAR && __column_numbers[LSTAR] != -1)
                    ||
                    (i == LSTAR && __column_numbers[LOG_LSTAR] != -1)
                )
            ) {
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
        std::valarray< std::list<double> > track_columns = parse_columns(
            track,
            header.get_all_columns()
        );
        
#ifndef NDEBUG
        track_columns[MTRACK].unique();
        assert(track_columns.size() == 1);
#endif
        __mass_list.push_back(track_columns[MTRACK].front());
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
                         double smooth_conv_inertia,
                         double smooth_rad_inertia,
                         double smooth_rad_mass,
                         int conv_inertia_nodes,
                         int rad_inertia_nodes,
                         int rad_mass_nodes) :
        track_quantities(NUM_COLUMNS)
    {
        std::cout << "Reading tracks from " << std::endl;
        DIR *dirstream=opendir(model_directory.c_str());
        std::string join;
        if(model_directory[model_directory.size()-1]=='/') join="";
        else join="/";
        if(dirstream==NULL)
            throw Error::PathNotFound("in MESA::Evolution constructor.",
                                      model_directory);
        struct dirent *entry; 
        while((entry=readdir(dirstream))) {
            std::string fname(entry->d_name);
            if(fname[0]!='.' && fname.substr(fname.size()-5)==".data") {
                std::cout << "Reading " << model_directory+join+entry->d_name
                    << std::endl;
                read_model_file(model_directory+join+entry->d_name);
            }
        }
        if(closedir(dirstream)) throw Error::Runtime(
                "Failed to close directory stream tied to "+model_directory+
                " in MESA::Evolution constructor.");
        sort_masses();
        std::valarray<double> masses = list_to_valarray(__mass_list);
        std::cout << "Done reading tracks." << std::endl
            << "Starting interpolation" << std::endl;
        interpolate_from(masses,
                         track_quantities[AGE],
                         track_quantities[RSTAR],
                         track_quantities[ICONV],
                         track_quantities[IRAD],
                         track_quantities[MRAD],
                         track_quantities[RRAD],
                         NaN,
                         smooth_conv_inertia,
                         smooth_rad_inertia,
                         smooth_rad_mass,
                         NaN,
                         0,
                         conv_inertia_nodes,
                         rad_inertia_nodes,
                         rad_mass_nodes,
                         0,
                         track_quantities[LSTAR],
                         NaN,
                         0);
        std::cout << "Done with interpolation" << std::endl;
    }
};
