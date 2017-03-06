#include "../StellarEvolution/CInterface.h"
#include "../Core/Common.h"

int main(int, char **)
{
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
