#include "MESAIO.h"

std::ostream &operator<<(std::ostream &os,
                         StellarEvolution::MESA::Column col)
{
    switch(col) {
        case StellarEvolution::MESA::AGE : os << "AGE"; break;
        case StellarEvolution::MESA::LOG_RSTAR : os << "LOG_RSTAR"; break;
        case StellarEvolution::MESA::LOG_LSTAR : os << "LOG_LSTAR"; break;
        case StellarEvolution::MESA::MRAD : os << "MRAD"; break;
        case StellarEvolution::MESA::RRAD : os << "RRAD"; break;
        case StellarEvolution::MESA::ICONV : os << "ICONV"; break;
        case StellarEvolution::MESA::IRAD : os << "IRAD"; break;
        default : throw Core::Error::IO("Unrecognized MESA column "
                                        "encountered!");
    };
    return os;
}

namespace StellarEvolution {
    namespace MESA {

        ///\brief A scaling constant used when transforming between different
        ///metallicity quantities.
        const double scaling = Yprotosun - Yprimordial + Zprotosun;

        double metallicity_from_feh(double feh)
        {
            return (
                feh
                +
                std::log10(
                    (1 - Yprimordial)
                    /
                    (Xprotosun + scaling * std::pow(10.0, feh))
                )
            );
        }

        double feh_from_metallicity(double metallicity)
        {
            return (
                metallicity
                -
                std::log10(
                    (1.0 - Yprimordial - scaling * std::pow(10.0, metallicity))
                    /
                    Xprotosun
                )
            );
        }

        const std::vector<QuantityID> Interpolator::__column_to_quantity(
                {
                StellarEvolution::NUM_QUANTITIES, //MTRACK - no quantity
                StellarEvolution::NUM_QUANTITIES, //AGE - no quantity
                StellarEvolution::RADIUS,         //LOG_RSTAR - converted
                StellarEvolution::RADIUS,         //RSTAR
                StellarEvolution::LUM,            //LOG_LSTAR - converted
                StellarEvolution::LUM,            //LSTAR
                StellarEvolution::MRAD,           //MRAD
                StellarEvolution::RRAD,           //RRAD
                StellarEvolution::ICONV,          //ICONV
                StellarEvolution::IRAD            //IRAD
                }
            );

        const std::vector<double> Interpolator::__default_smoothing(
            {
            Core::NaN,  //RADIUS
            5,          //ICONV
            Core::NaN,  //LUM
            6,          //IRAD
            7,          //MRAD
            6           //RRAD: TBD
            }
        );

        const std::vector<int> Interpolator::__default_nodes(
            {
            0,      //RADIUS
            3000,   //ICONV
            0,      //LUM
            3000,   //IRAD
            6000,   //MRAD
            3000    //RRAD: TBD
            }
        );

        const std::vector<bool> Interpolator::__default_vs_log_age(
            {
            true,   //RADIUS
            true,   //ICONV
            true,   //LUM
            true,   //IRAD
            true,   //MRAD
            true    //RRAD
            }
        );

        const std::vector<bool> Interpolator::__default_log_quantity(
            {
            false,  //RADIUS
            false,   //ICONV
            false,  //LUM
            false,  //IRAD
            false,  //MRAD
            false   //RRAD
            }
        );

        void Header::set_column_names()
        {
            __column_names.resize(NUM_COLUMNS);
            __column_names[MESA::MTRACK] = "M_ini";
            __column_names[MESA::AGE] = "age";
            __column_names[MESA::LOG_RSTAR] = "";
            __column_names[MESA::RSTAR] = "R_star";
            __column_names[MESA::LOG_LSTAR] = "";
            __column_names[MESA::LSTAR] = "L_star";
            __column_names[MESA::MRAD] = "M_rad";
            __column_names[MESA::RRAD] = "R_tachocline";
            __column_names[MESA::ICONV] = "I_conv";
            __column_names[MESA::IRAD] = "I_rad";
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
                    throw Core::Error::IO(msg.str());
                }
                std::getline(track, line);
            }

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
            __column_numbers.resize(MESA::NUM_COLUMNS, -1);

            unsigned line_number=0;
            read_column_numbers(track, filename, line_number);

            for(int i=0; i<MESA::NUM_COLUMNS; ++i)
                if(
                    __column_numbers[i]==-1
                    &&
                    !(
                        (
                            i == MESA::LOG_RSTAR
                            &&
                            __column_numbers[MESA::RSTAR] != -1
                        ) 
                        ||
                        (
                            i == MESA::RSTAR
                            &&
                            __column_numbers[MESA::LOG_RSTAR] != -1
                        )
                        ||
                        (
                            i == MESA::LOG_LSTAR
                            &&
                            __column_numbers[MESA::LSTAR] != -1
                        )
                        ||
                        (
                            i == MESA::LSTAR
                            &&
                            __column_numbers[MESA::LOG_LSTAR] != -1
                        )
                    )
                ) {
                    std::ostringstream message;
                    message << "Failed to find '" << __column_names[i]
                        << "' on line " << line_number << " of " << filename
                        << ".";
                    throw Core::Error::IO(message.str());
                }
        }

        EvolutionIterator &EvolutionIterator::operator=(
            const EvolutionIterator &rhs)
        {
            mass_iter = rhs.mass_iter;
            feh_iter = rhs.feh_iter;
            age_iter = rhs.age_iter;
            quantity_iter = rhs.quantity_iter;
            return *this;
        }

        EvolutionIterator &EvolutionIterator::operator++()
        {
            ++mass_iter;
            ++feh_iter;
            ++age_iter;
            for(size_t i=0; i<quantity_iter.size(); ++i) ++quantity_iter[i];
            return *this;
        }

#ifndef NDEBUG
        void Interpolator::log_current_age_ranges() const
        {
            assert(__mass_list.size() == __feh_list.size());
            std::list<double>::const_iterator 
                mass_iter = __mass_list.begin(),
                feh_iter = __feh_list.begin();

            std::list< std::valarray<double> >::const_iterator
                age_iter = __track_ages.begin();
            while(mass_iter != __mass_list.end()) {
                assert(feh_iter != __feh_list.end());
                assert(age_iter != __track_ages.end());
                ++mass_iter;
                ++feh_iter;
                ++age_iter;
            }
        }
#endif

        bool Interpolator::parse_model_file_name(const std::string &filename)
        {
            const double PRECISION = 1000.0;
            std::istringstream fname_stream(filename);
            if(fname_stream.get() != 'M' || !fname_stream.good())
                return false;

            double mass = 0;
            fname_stream >> mass;
            if(!fname_stream.good() || mass <= 0) return false;

            if(fname_stream.get() != '_' || !fname_stream.good())
                return false;
            if(fname_stream.get() != 'Z' || !fname_stream.good())
                return false;

            double metallicity = 0;
            fname_stream >> metallicity;
            if(!fname_stream.good() || metallicity <= 0) return false;

            __feh_list.push_back(
                std::round(
                    PRECISION
                    *
                    feh_from_metallicity(std::log10(metallicity / Zprotosun))
                )
                /
                PRECISION
            );
            __mass_list.push_back(std::round(PRECISION * mass) / PRECISION);
            return true;
        }

        ///Used as comparison when sorting quantities by age.
        class CompareAges {
        private:
            const std::valarray<double> __ages;
        public:
            CompareAges(const std::valarray<double> &ages) : __ages(ages) {}

            bool operator()(size_t i1, size_t i2)
            {return __ages[i1] < __ages[i2];}
        };

        void Interpolator::sort_last_track_by_age()
        {
            if(
                std::is_sorted(std::begin(__track_ages.back()), 
                               std::end(__track_ages.back()))
            )
                return;

            std::valarray<size_t> age_sorting_indices(
                __track_ages.back().size()
            );
            for(size_t i = 0; i < __track_ages.back().size(); ++i)
                age_sorting_indices[i] = i;
            std::sort(std::begin(age_sorting_indices),
                      std::end(age_sorting_indices),
                      CompareAges(__track_ages.back()));

            __track_ages.back() = std::valarray<double>(
                __track_ages.back()[age_sorting_indices]
            );
            assert(
                std::is_sorted(std::begin(__track_ages.back()),
                               std::end(__track_ages.back()))
            );
            for(
                int quantity = 0;
                quantity < StellarEvolution::NUM_QUANTITIES;
                ++quantity
            ) 
                __track_quantities[quantity].back() = std::valarray<double>(
                    __track_quantities[quantity].back()[age_sorting_indices]
                );
        }

        void Interpolator::read_model_file(const std::string &filename)
        {
            std::ifstream track(filename.c_str());
            Header header(track, filename);
            std::vector< std::list<double> > track_columns = parse_columns(
                track,
                header.get_all_columns()
            );

#ifndef NDEBUG
            track_columns[MESA::MTRACK].unique();
            assert(track_columns[MESA::MTRACK].size() == 1);
            assert(__mass_list.back()
                   ==
                   track_columns[MESA::MTRACK].front());
#endif
            __track_ages.push_back(
                Core::list_to_valarray(track_columns[MESA::AGE]) * 1e-9
            );
            std::clog
                << "Model: M = " << __mass_list.back()
                << ", [Fe/H] = " << __feh_list.back()
                << ", t_max = "
                << __track_ages.back()[__track_ages.back().size() - 1]
                << std::endl;
            for(int column = 0; column < MESA::NUM_COLUMNS; ++column) {
                QuantityID quantity = __column_to_quantity[column];
                if(
                    quantity == StellarEvolution::NUM_QUANTITIES
                    ||
                    track_columns[column].empty()
                ) continue;

                if(column == MESA::LOG_RSTAR || column == MESA::LOG_LSTAR)
                    __track_quantities[quantity].push_back(
                        std::pow(
                            10.0,
                            Core::list_to_valarray(track_columns[column])
                        )
                    );
                else 
                    __track_quantities[quantity].push_back(
                        Core::list_to_valarray(track_columns[column])
                    );
            }
            sort_last_track_by_age();
        }

        void Interpolator::get_mass_feh_grid(
            std::valarray<double> &masses,
            std::valarray<double> &feh
        )
        {
            assert(__mass_list.size() == __feh_list.size());
            std::list<double>::const_iterator 
                mass_iter = __mass_list.begin(),
                feh_iter = __feh_list.begin();

            std::list< std::valarray<double> >::const_iterator
                age_iter = __track_ages.begin();

            size_t num_masses = 0;
            for(
                double first_feh = *feh_iter;
                *feh_iter == first_feh;
                ++feh_iter
            )
                ++num_masses;
            feh_iter = __feh_list.begin();
            assert(__mass_list.size() % num_masses == 0);
            size_t num_feh = __mass_list.size() / num_masses;

            masses.resize(num_masses);
            feh.resize(num_feh);

            for(size_t feh_index = 0; feh_index < num_feh; ++feh_index) {
                for(
                    size_t mass_index = 0;
                    mass_index < num_masses;
                    ++mass_index
                ) {
                    assert(mass_iter != __mass_list.end());
                    assert(feh_iter != __feh_list.end());
                    std::cerr
                        << "M = " << *mass_iter
                        << ", [Fe/H] = " << *feh_iter
                        << ", age range = " << (*age_iter)[0]
                        << " - " << (*age_iter)[age_iter->size() -1]
                        << std::endl;
                    if(feh_index == 0) 
                        masses[mass_index] = *mass_iter;
                    else if(masses[mass_index] != *mass_iter)
                        throw Core::Error::IO(
                            "Input MESA tracks do not lie on a grid of "
                            "masses and [Fe/H]."
                        );

                    if(mass_index == 0) 
                        feh[feh_index] = *feh_iter;
                    else if(feh[feh_index] != *feh_iter)
                        throw Core::Error::IO(
                            "Input MESA tracks do not lie on a grid of "
                            "masses and [Fe/H]."
                        );

                    ++mass_iter;
                    ++feh_iter;
                    ++age_iter;
                }
            }
            assert(feh_iter == __feh_list.end());
            assert(mass_iter == __mass_list.end());
        }

        EvolutionIterator Interpolator::begin()
        {
            EvolutionIterator result;
            result.mass_iter = __mass_list.begin();
            result.feh_iter = __feh_list.begin();
            result.age_iter = __track_ages.begin();
            for(size_t quantity = 0; quantity < NUM_QUANTITIES; ++quantity)
                result.quantity_iter[quantity] =
                    __track_quantities[quantity].begin();
            return result;
        }

        EvolutionIterator Interpolator::end()
        {
            EvolutionIterator result;
            result.mass_iter = __mass_list.end();
            result.feh_iter = __feh_list.end();
            result.age_iter = __track_ages.end();
            for(size_t quantity = 0; quantity < NUM_QUANTITIES; ++quantity)
                result.quantity_iter[quantity] =
                    __track_quantities[quantity].end();
            return result;
        }

        void Interpolator::move(EvolutionIterator &dest,
                                EvolutionIterator &source)
        {
            __mass_list.splice(dest.mass_iter,
                               __mass_list,
                               source.mass_iter);

            __feh_list.splice(dest.feh_iter,
                              __feh_list,
                              source.feh_iter);

            __track_ages.splice(dest.age_iter,
                                __track_ages,
                                source.age_iter);

            for(size_t quantity = 0; quantity < NUM_QUANTITIES; ++quantity)
                __track_quantities[quantity].splice(
                    dest.quantity_iter[quantity],
                    __track_quantities[quantity],
                    source.quantity_iter[quantity]
                );
#ifndef NDEBUG
            log_current_age_ranges();
#endif
        }

        void Interpolator::sort_tracks()
        {
            EvolutionIterator iter = begin(), stop_iter = end();
            double last_sorted_mass = *iter.mass_iter,
                   last_sorted_feh = *iter.feh_iter;
            while(iter != stop_iter) {
                while(
                    iter != stop_iter
                    &&
                    (
                        *iter.feh_iter > last_sorted_feh
                        ||
                        (
                            *iter.feh_iter == last_sorted_feh
                            &&
                            *iter.mass_iter >= last_sorted_mass
                        )
                    )
                ) {
                    last_sorted_mass = *iter.mass_iter;
                    last_sorted_feh = *iter.feh_iter;
                    ++iter;
                }
                if(iter == stop_iter) break;
                EvolutionIterator dest = begin();

                while(
                    *(dest.feh_iter) < *(iter.feh_iter)
                    ||
                    (
                        *(dest.feh_iter) == *(iter.feh_iter)
                        &&
                        *(dest.mass_iter) <= *(iter.mass_iter)
                    )
                )
                    ++dest;
                EvolutionIterator source = iter++;
                move(dest, source);
            }
        }

        Interpolator::Interpolator(const std::string &model_directory,
                                   unsigned num_threads,
                                   const std::vector<double> &smoothing,
                                   const std::vector<int> &nodes,
                                   const std::vector<bool> &vs_log_age,
                                   const std::vector<bool> &log_quantity) :
            __track_quantities(NUM_QUANTITIES)
        {
            std::clog << "Reading tracks from "
                << model_directory
                << std::endl;
            DIR *dirstream = opendir(model_directory.c_str());
            std::string join;
            if(model_directory[model_directory.size()-1] == '/') join="";
            else join = "/";
            if(dirstream == NULL)
                throw Core::Error::PathNotFound(
                    "in MESA::Evolution constructor.",
                    model_directory
                );
            struct dirent *entry;
            while((entry = readdir(dirstream))) {
                std::string fname(entry->d_name);
                if(
                    fname[0] != '.'
                    &&
                    fname.substr(fname.size()-4) == ".csv"
                    &&
                    parse_model_file_name(fname)
                ) {
                    std::cout
                        << "Reading "
                        << model_directory + join + fname
                        << std::endl;
                    read_model_file(model_directory + join + fname);
                } else
                    std::cout
                        << "Skipping "
                        << model_directory + join + entry->d_name
                        << std::endl;
            }
            if(closedir(dirstream)) throw Core::Error::Runtime(
                "Failed to close directory stream tied to "
                +
                model_directory
                +
                " in MESA::Evolution constructor."
            );
            sort_tracks();
            std::valarray<double> masses, feh;
            get_mass_feh_grid(masses, feh);
            std::cout
                << "Done reading tracks." << std::endl
                << "Starting interpolation" << std::endl;

            create_from(masses,
                        feh,
                        __track_ages,
                        __track_quantities,
                        smoothing,
                        nodes,
                        vs_log_age,
                        log_quantity,
                        num_threads);
            std::cout << "Done with interpolation" << std::endl;
        }
    } //End MESA namespace.
} // End StellarEvolution namespace.
