#include "test_MESAIO.h"

test_MESAIO::test_MESAIO(const std::string& mesa_dir) :
    __mesa_dir(mesa_dir) 
{
    TEST_ADD(test_MESAIO::test_single_tracks);
}

void test_MESAIO::test_single_tracks()
{
    std::cout << "Starting test" << std::endl;
    const std::string& filename = "../MESA/FeH=+0.0/M1.0_Z0.015.csv";
    std::ifstream track(filename.c_str());
    MESA::Header header(track, filename);
    std::valarray< std::list<double> > track_columns = parse_columns(
        track,
        header.get_all_columns()
    );
    MESA::Evolution evolution("../MESA/FeH=+0.0");
    const EvolvingStellarQuantity
        *R_star = evolution.interpolate_radius(1.0),
        *L_star = evolution.interpolate_luminosity(1.0),
        *M_core = evolution.interpolate_zone_mass(1.0, radiative),
        *R_core = evolution.interpolate_core_boundary(1.0),
        *Iconv_star = evolution.interpolate_moment_of_inertia(1.0, convective),
        *Irad_star = evolution.interpolate_moment_of_inertia(1.0, radiative);

    std::cout << "Age" << std::endl << std::endl;
    typedef std::list<double>::const_iterator TrackIter;
    for(
        TrackIter age_i = track_columns[MESA::AGE].begin();
        age_i != track_columns[MESA::AGE].end();
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
            std::cout
                << age << "," 
                << (*R_star)(age) << ","
                << (*L_star)(age) << ","
                << 1.0 << ","
                << (*R_core)(age) << ","
                << 1.0 << ","
                << 1.0 - (*M_core)(age) << ","
                << (*M_core)(age) << ","
                << (*Iconv_star)(age) << ","
                << (*Irad_star)(age)
                << std::endl;
        }
    }
    delete R_star;
    delete L_star;
    delete M_core;
    delete R_core;
    delete Iconv_star;
    delete Irad_star;
}

#ifdef STANDALONE
int main()
{
    std::srand(std::time(NULL));
    std::cout.setf(std::ios_base::scientific);
    std::cout.precision(16);
    Test::TextOutput output(Test::TextOutput::Verbose);
    test_MESAIO tests;
    return (tests.run(output, false) ? EXIT_SUCCESS : EXIT_FAILURE);
}
#endif
