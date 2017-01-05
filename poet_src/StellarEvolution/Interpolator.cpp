/**\file
 *
 * \brief Defines some of the methods of the StellarEvolution class.
 * 
 * \ingroup StellarSystem_group
 */

#include "Interpolator.h"
#include "SumQuantity.h"
#include "../Core/Functions.h"

namespace StellarEvolution {

    int Interpolator::find_first_core_index(
        const std::valarray<double> &core_mass
    ) const
    {
        int first_core_index = 0;
        while(
            first_core_index < static_cast<int>(core_mass.size())
            &&
            core_mass[first_core_index] == 0
        ) ++first_core_index;
        return first_core_index;
    }

    void Interpolator::create_from(
        const std::valarray<double> &tabulated_masses,
        const std::valarray<double> &tabulated_metallicities,
        const std::list< std::valarray<double> > &tabulated_ages,
        const std::vector< std::list< std::valarray<double> > >
        &tabulated_quantities,
        const std::vector<double> &smoothing,
        const std::vector<int> &nodes,
        const std::vector<bool> &vs_log_age,
        const std::vector<bool> &log_quantity
    )
    {
        __track_masses = tabulated_masses;
        __track_metallicities = tabulated_metallicities;
        __core_formation = Core::Inf;
        __vs_log_age = vs_log_age;
        __log_quantity = log_quantity;

        size_t num_tracks = (tabulated_masses.size()
                             *
                             tabulated_metallicities.size());

        assert(tabulated_ages.size() == num_tracks);
        assert(tabulated_quantities.size() == NUM_QUANTITIES);
        assert(nodes.size() == NUM_QUANTITIES);
        assert(smoothing.size() == NUM_QUANTITIES);
        assert(vs_log_age.size() == NUM_QUANTITIES);
        assert(log_quantity.size() == NUM_QUANTITIES);

        typedef std::list< std::valarray<double> >::const_iterator 
            track_quantity_iter;

        track_quantity_iter ages_iter=tabulated_ages.begin();

        std::vector< track_quantity_iter > track_iter(NUM_QUANTITIES);

        for(
            int quantity_index = 0;
            quantity_index < NUM_QUANTITIES;
            ++quantity_index
        ) {
            assert(tabulated_quantities[quantity_index].size()
                   ==
                   num_tracks);
            __interpolated_quantities[quantity_index].resize(num_tracks);
            track_iter[quantity_index] = 
                tabulated_quantities[quantity_index].begin();
        }

        for(size_t grid_index = 0; grid_index < num_tracks; ++grid_index) {

            std::clog << "Grid point " << grid_index
                << ": M = "
                << tabulated_masses[grid_index % tabulated_masses.size()]
                << ", [Fe/H] = " 
                << tabulated_metallicities[grid_index
                                           /
                                           tabulated_masses.size()]
                << std::endl;

            int first_core_index = (
                tabulated_quantities[MRAD].empty()
                ? 0
                : find_first_core_index(*(track_iter[MRAD]))
            );

            bool no_core = (
                first_core_index == static_cast<int>(
                    track_iter[MRAD]->size()
                )
                ||
                (*(track_iter[MRAD]))[first_core_index] == 0
            );

            if(!no_core) 
                __core_formation = std::min(
                    __core_formation,
                    (*ages_iter)[std::max(0, first_core_index - 1)]
                );

            std::valarray<double> log_ages = std::log(*ages_iter);

            for(
                int quantity_index = 0;
                quantity_index < NUM_QUANTITIES;
                ++quantity_index
            ) {
                std::clog << "Interpolating "
                    << static_cast<QuantityID>(quantity_index)
                    << std::endl;

                if(quantity_index >= FIRST_CORE_QUANTITY && no_core) {
                    __interpolated_quantities[quantity_index][grid_index] = 
                        new Core::ZeroFunction();
                    __vs_log_age[quantity_index] = false;
                    __log_quantity[quantity_index] = false;
                    continue;
                }

                int first_interp_index;
                const std::valarray<double> *interp_quantity;
                if(log_quantity[quantity_index])
                    interp_quantity = new std::valarray<double>(
                        std::log(*(track_iter[quantity_index]))
                    );
                else
                    interp_quantity = &(*(track_iter[quantity_index]));
                if(quantity_index < FIRST_CORE_QUANTITY)
                    first_interp_index = 0;
                else if(log_quantity[quantity_index])
                    first_interp_index = first_core_index;
                else
                    first_interp_index = std::max(0, first_core_index - 1);

                __interpolated_quantities[quantity_index][grid_index] = 
                    new Core::InterpolatingFunctionALGLIB(
                        (
                            vs_log_age[quantity_index]
                            ? &(log_ages[first_interp_index])
                            : &((*ages_iter)[first_interp_index])
                        ),
                        &((*interp_quantity)[first_interp_index]),
                        log_ages.size() - first_interp_index,
                        NULL,
                        smoothing[quantity_index],
                        nodes[quantity_index]
                    );
                if(log_quantity[quantity_index]) delete interp_quantity;
            }

            ++ages_iter;
            for(
                size_t quantity_index = 0; 
                quantity_index < NUM_QUANTITIES;
                ++quantity_index
            ) ++track_iter[quantity_index];
        }
    }

    EvolvingStellarQuantity *Interpolator::operator()(
        QuantityID quantity,
        double mass,
        double metallicity
    ) const
    {
        return new EvolvingStellarQuantity(
            mass,
            metallicity,
            __track_masses, 
            __track_metallicities,
            __interpolated_quantities[quantity],
            __vs_log_age[quantity],
            __log_quantity[quantity],
            quantity >= FIRST_CORE_QUANTITY
        );
    }

#ifndef NO_SERIALIZE
    void Interpolator::load_state(const std::string &filename)
    {
        std::ifstream ifs(filename.c_str());
        boost::archive::text_iarchive ia(ifs);
        ia >> (*this);
        ifs.close();
    }

    void Interpolator::save_state(const std::string &filename) const
    {
        std::ofstream ofs(filename.c_str());
        boost::archive::text_oarchive oa(ofs);
        oa << (*this);
    }

} //End StellarEvolution namespace.
#endif
