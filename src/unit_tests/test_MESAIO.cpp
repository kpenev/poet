#include "test_MESAIO.h"

test_MESAIO::test_MESAIO(const std::string& mesa_dir) :
    __mesa_dir(mesa_dir) 
{
    TEST_ADD(test_MESAIO::test_single_tracks);
}

void test_MESAIO::test_single_tracks()
{
    std::cout << "Starting test" << std::endl;
    const std::string& filename = "../MESA/0.7Msun_table.csv";
    std::ifstream track(filename.c_str());
    MESA::Header header(track, filename);
    std::valarray< std::list<double> > track_columns = parse_columns(
        track,
        header.get_all_columns()
    );
    TEST_THROWS_NOTHING(MESA::Evolution evolution("../MESA"));

    std::cout << "Age" << std::endl << std::endl;
    typedef std::list<double>::const_iterator TrackIter;
    for(
        TrackIter age_i = track_columns[MESA::AGE].begin();
        age_i != track_columns[MESA::AGE].end();
        ++age_i
    ) std::cout << *age_i << std::endl;
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
