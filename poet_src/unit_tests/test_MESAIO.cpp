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
    Interpolator evolution(
        "../StellarEvolution/MESA"
    );
    const std::vector<std::string> quantity_name({"R*",
                                                  "L*",
                                                  "Mcore",
                                                  "Rcore",
                                                  "Ienv",
                                                  "Icore"});
    const std::vector<const EvolvingStellarQuantity*> quantity(
        {
        evolution.interpolate_radius(1.0, 0.0),
        evolution.interpolate_luminosity(1.0, 0.0),
        evolution.interpolate_core_mass(1.0, 0.0),
        evolution.interpolate_core_boundary(1.0, 0.0),
        evolution.interpolate_moment_of_inertia(1.0, 0.0, convective),
        evolution.interpolate_moment_of_inertia(1.0, 0.0, radiative)
        }
    );

    double start_age = (track_columns[StellarEvolution::MESA::AGE].front() 
                        *
                        1e-9);
    std::cout << "Age";
    for(
        size_t quantity_ind = 0;
        quantity_ind < quantity.size();
        ++quantity_ind
    ) {
        std::cout << "," << quantity_name[quantity_ind];
        quantity[quantity_ind]->select_interpolation_region(start_age);
    }
    std::cout << std::endl;

    typedef std::list<double>::const_iterator TrackIter;
    for(
        TrackIter age_i = track_columns[StellarEvolution::MESA::AGE].begin();
        age_i != track_columns[StellarEvolution::MESA::AGE].end();
        ++age_i
    ) {
        TrackIter next_age = age_i;
        ++next_age;
        for(int substep = 0; substep < 10; ++substep) {
            double age = (
                (*age_i) * (10 - substep)
                + 
                (*next_age) * substep
            ) / 1e10;
            std::cout << age << ",";
            for(
                size_t quantity_ind = 0;
                quantity_ind < quantity.size();
                ++quantity_ind
            ) {

                if(age > quantity[quantity_ind]->next_discontinuity())
                    quantity[quantity_ind]->enable_next_interpolation_region();
                std::cout
                    << (*(quantity[quantity_ind]))(age) << ",";
            }
            std::cout << std::endl;
        }
    }

    for(
        size_t quantity_ind = 0;
        quantity_ind < quantity.size();
        ++quantity_ind
    )
        delete quantity[quantity_ind];
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
