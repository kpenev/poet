#include "test_MESAIO.h"

test_MESAIO::test_MESAIO(const std::string& mesa_dir) :
    __mesa_dir(mesa_dir) 
{
    TEST_ADD(test_MESAIO::test_single_tracks);
}

void test_MESAIO::test_single_tracks()
{
    using StellarEvolution::MESA::Interpolator;
    using StellarEvolution::EvolvingStellarQuantity;
    using Core::convective;
    using Core::radiative;

    std::cout << "Starting test" << std::endl;
    const std::string& filename = "../StellarEvolution/MESA/M1.0_Z0.015.csv";
    std::ifstream track(filename.c_str());
    StellarEvolution::MESA::Header header(track, filename);
    std::vector< std::list<double> > track_columns = parse_columns(
        track,
        header.get_all_columns()
    );
    track_columns[StellarEvolution::MESA::AGE].sort();

    Interpolator saved_evolution(
        "../../MESA_tracks",
        {0.1, 0.2, 0.3, 0.4, 0.5, 0.6}, //smoothing
        {100, 100, 100, 100, 100, 100}, //nodes
        {true, true, true, true, true, true}, //vs log(age)
        {true, false, true, false, true, false} //log(quantity)
    );
    saved_evolution.save_state("serialized_evolution");
    Interpolator loaded_evolution;
    loaded_evolution.load_state("serialized_evolution");
/*    evolution.load_state(
        "../../stellar_evolution_interpolators/"
        "0d728daf-63a2-419c-b254-51e8f1c82803"
    );*/
    const std::vector<std::string> quantity_name({"R*",
                                                  "L*",
                                                  "Mcore",
                                                  "Rcore",
                                                  "Ienv",
                                                  "Icore"});
    const std::vector<const EvolvingStellarQuantity*> saved_quantity(
        {
        saved_evolution(StellarEvolution::RADIUS, 1.0 + 1e-5, 0.0 + 1e-5),
        saved_evolution(StellarEvolution::LUM, 1.0 + 1e-5, 0.0 + 1e-5),
        saved_evolution(StellarEvolution::MRAD, 1.0 + 1e-5, 0.0 + 1e-5),
        saved_evolution(StellarEvolution::RRAD, 1.0 + 1e-5, 0.0 + 1e-5),
        saved_evolution(StellarEvolution::ICONV, 1.0 + 1e-5, 0.0 + 1e-5),
        saved_evolution(StellarEvolution::IRAD, 1.0 + 1e-5, 0.0 + 1e-5)
        }
    );
    const std::vector<const EvolvingStellarQuantity*> loaded_quantity(
        {
        loaded_evolution(StellarEvolution::RADIUS, 1.0 + 1e-5, 0.0 + 1e-5),
        loaded_evolution(StellarEvolution::LUM, 1.0 + 1e-5, 0.0 + 1e-5),
        loaded_evolution(StellarEvolution::MRAD, 1.0 + 1e-5, 0.0 + 1e-5),
        loaded_evolution(StellarEvolution::RRAD, 1.0 + 1e-5, 0.0 + 1e-5),
        loaded_evolution(StellarEvolution::ICONV, 1.0 + 1e-5, 0.0 + 1e-5),
        loaded_evolution(StellarEvolution::IRAD, 1.0 + 1e-5, 0.0 + 1e-5)
        }
    );

    double start_age = (track_columns[StellarEvolution::MESA::AGE].front() 
                        *
                        1e-11);
    std::cout << "Age";
    for(
        size_t quantity_ind = 0;
        quantity_ind < saved_quantity.size();
        ++quantity_ind
    ) {
        std::cout << "," << quantity_name[quantity_ind];
        saved_quantity[quantity_ind]->select_interpolation_region(
            std::max(start_age, saved_quantity[quantity_ind]->range_low())
        );
        loaded_quantity[quantity_ind]->select_interpolation_region(
            std::max(start_age, saved_quantity[quantity_ind]->range_low())
        );
    }
    std::cout << std::endl;

    track_columns[StellarEvolution::MESA::AGE].push_front(start_age * 1e9);
    typedef std::list<double>::const_iterator TrackIter;
    TrackIter next_age = track_columns[StellarEvolution::MESA::AGE].begin();
    ++next_age;
    for(
        TrackIter age_i = track_columns[StellarEvolution::MESA::AGE].begin();
        next_age != track_columns[StellarEvolution::MESA::AGE].end();
        ++age_i
    ) {
        assert(*next_age > *age_i);
        for(int substep = 0; substep < 10; ++substep) {
            double age = (
                (*age_i) * (10 - substep)
                + 
                (*next_age) * substep
            ) / 1e10;
            std::cout << age << ",";
            for(
                size_t quantity_ind = 0;
                quantity_ind < saved_quantity.size();
                ++quantity_ind
            ) {
                if(age < saved_quantity[quantity_ind]->range_low()) {
                    std::cout << Core::NaN << ",";
                    continue;
                }
                if(age > saved_quantity[quantity_ind]->range_high()) {
                    std::cout << "*,";
                    continue;
                }

                while(
                    age
                    >
                    saved_quantity[quantity_ind]->next_discontinuity()
                ) {
                    assert(
                        age
                        >
                        loaded_quantity[quantity_ind]->next_discontinuity()
                    );
                    (
                        saved_quantity[quantity_ind]
                    )->enable_next_interpolation_region();
                    (
                        loaded_quantity[quantity_ind]
                    )->enable_next_interpolation_region();
                }
                const Core::FunctionDerivatives
                    *deriv = saved_quantity[quantity_ind]->deriv(age),
                    *loaded_deriv = loaded_quantity[quantity_ind]->deriv(age);
                double value = deriv->order(0);
                assert(value == loaded_deriv->order(0));
                std::cout << value << ",";
                delete deriv;
                delete loaded_deriv;
            }
            std::cout << std::endl;
        }
        ++next_age;
    }

    for(
        size_t quantity_ind = 0;
        quantity_ind < saved_quantity.size();
        ++quantity_ind
    ) {
        delete saved_quantity[quantity_ind];
        delete loaded_quantity[quantity_ind];
    }
    saved_evolution.delete_tracks();
    loaded_evolution.delete_tracks();
}

#ifdef STANDALONE
int main()
{
    std::srand(std::time(NULL));
    std::cout.setf(std::ios_base::scientific);
    std::cout.precision(16);
    Test::TextOutput output(Test::TextOutput::Verbose);
    test_MESAIO tests;
    tests.test_single_tracks();
    return (tests.run(output, false) ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
