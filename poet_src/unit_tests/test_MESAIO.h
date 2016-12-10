/**\file 
 *
 * \brief Declares the test suite for exercising the MESAIO class.
 *
 * \ingroup UnitTests_group
 */

#ifndef __TEST_MESAIO_H
#define __TEST_MESAIO_H

//#include "Common.h"
#include "StellarEvolution/interface/MESAIO.h"
#include "IO/IOUtil.h"
#include <cpptest.h>
#include <sstream>
#include <cmath>
#include <ctime>

class test_MESAIO : public Test::Suite {
private:
    ///The directory where the MESA tracks are located.
    std::string __mesa_dir;

protected:
    ///No fixtures at this time
    void setup() {};

    ///No fixtures at this time
    void tear_down() {};

public:
    ///Test the interpolation of stellar evolution based on MESA tracks.
    test_MESAIO(
        ///The directory from which to read MESA files.
        const std::string& mesa_dir = "../StellarEvolution/MESA/"
    );

    ///Test a single track's interpolation.
    void test_single_tracks();
};

#endif
