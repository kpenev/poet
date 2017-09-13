#include "../StellarEvolution/CInterface.h"
#include "../Core/Common.h"

int main(int, char **)
{
    std::cout << "[Z](feh = -1.01444556) = "
              << metallicity_from_feh(-1.01444556)
              << std::endl;
    std::cout << "[Z](feh = -0.76321518) = "
              << metallicity_from_feh(-0.76321518)
              << std::endl;
    std::cout << "[Z](feh = -0.51101857) = "
              << metallicity_from_feh(-0.51101857)
              << std::endl;
    std::cout << "[Z](feh = -0.25708473) = "
              << metallicity_from_feh(-0.25708473)
              << std::endl;
    std::cout << "[Z](feh = 0) = "
              << metallicity_from_feh(0)
              << std::endl;
    std::cout << "[Z](feh = 0.2628914) = "
              << metallicity_from_feh(0.2628914)
              << std::endl;
    std::cout << "[Z](feh = 0.53680605) = "
              << metallicity_from_feh(0.53680605)
              << std::endl;

    for(double z = -1.0; z <= 0.6; z+=0.25)
        std::cout << "[Fe/H]([Z] = " << z << ") = "
                  << feh_from_metallicity(z)
                  << std::endl;

    double smoothing[] = {Core::NaN,
                          Core::NaN,
                          Core::NaN,
                          Core::NaN,
                          Core::NaN,
                          Core::NaN};
    int nodes[] = {0, 0, 0, 0, 0, 0};
    bool vs_log_age[] = {true, true, true, true, true, true};
    bool log_quantity[] = {false, false, false, false, false, false};
    MESAInterpolator *interpolator = create_interpolator(
        "../MESA_tracks",
        smoothing,
        nodes,
        vs_log_age,
        log_quantity,
        4
    );
    destroy_interpolator(interpolator);
}
