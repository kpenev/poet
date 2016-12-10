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

    size_t Interpolator::find_first_core_index(
        const std::valarray<double> &core_mass
    ) const
    {
        size_t first_core_index = 0;
        while(
            first_core_index < core_mass.size()
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
        const std::vector<int> &nodes
    )
    {
        __track_masses = tabulated_masses;
        __track_metallicities = tabulated_metallicities;
        __core_formation = Core::Inf;

        size_t num_tracks = (tabulated_masses.size()
                             *
                             tabulated_metallicities.size());
        assert(tabulated_ages.size() == num_tracks);
        assert(tabulated_quantities.size() == NUM_QUANTITIES);

        typedef 
            std::list< std::valarray<double> >::const_iterator 
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
            std::clog << "Track " << grid_index << std::endl;

            std::valarray<double> log_ages = std::log(*ages_iter);
            for(
                int quantity_index = 0;
                quantity_index < FIRST_CORE_QUANTITY;
                ++quantity_index
            ) {
                std::clog
                    << "Interpolating " << QUANTITY_NAME[quantity_index]
                    << "("
                    << "num_ages: " << log_ages.size()
                    << ", num_radii: " << track_iter[quantity_index]->size()
                    << ", smoothing: " << smoothing[quantity_index]
                    << ", nodes: " << nodes[quantity_index]
                    << ", range: " << log_ages[0] 
                    << " - " << log_ages[log_ages.size() - 1]
                    << ")...";
                std::clog.flush();

                __interpolated_quantities[quantity_index][grid_index] = 
                    new Core::InterpolatingFunctionALGLIB(
                        log_ages,
                        *(track_iter[quantity_index]),
                        std::valarray<double>(),
                        smoothing[quantity_index],
                        nodes[quantity_index]
                    );
                std::clog << "done" << std::endl;
            }

            size_t first_core_index = (
                tabulated_quantities[MRAD].empty()
                ? 0
                : find_first_core_index(*(track_iter[MRAD]))
            );
            if(
                first_core_index == track_iter[MRAD]->size()
                ||
                (*(track_iter[MRAD]))[first_core_index] == 0
            ) {
                for(
                    size_t quantity_index = FIRST_CORE_QUANTITY;
                    quantity_index < NUM_QUANTITIES;
                    ++quantity_index
                ) __interpolated_quantities[quantity_index][grid_index] = 
                    new Core::ZeroFunction();
            } else {
                if(first_core_index > 0) --first_core_index;
                std::slice core_slice(first_core_index,
                                      log_ages.size() - first_core_index,
                                      1);
                std::valarray<double> core_log_ages = log_ages[core_slice];
                __core_formation = std::min(__core_formation,
                                            (*ages_iter)[first_core_index]);
                std::clog 
                    << "core_formation = " << __core_formation 
                    << std::endl;
                for(
                    size_t quantity_index = FIRST_CORE_QUANTITY;
                    quantity_index < NUM_QUANTITIES;
                    ++quantity_index
                ) {
                    std::clog
                        << "Interpolating " 
                        << QUANTITY_NAME[quantity_index] << "("
                        << "num_ages: " << log_ages.size()
                        << ", num_points: " 
                        << track_iter[quantity_index]->size()
                        << ", smoothing: " << smoothing[quantity_index]
                        << ", nodes: " << nodes[quantity_index]
                        << ", range: " << core_log_ages[0] 
                        << " - " << core_log_ages[core_log_ages.size() - 1]
                        << ")...";
                    std::clog.flush();

                    __interpolated_quantities[quantity_index][grid_index] = 
                        new Core::InterpolatingFunctionALGLIB(
                            core_log_ages,
                            (*(track_iter[quantity_index]))[core_slice],
                            std::valarray<double>(),
                            smoothing[quantity_index],
                            nodes[quantity_index]
                        );
                    std::clog << "done" << std::endl;
                }
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
            true,
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
