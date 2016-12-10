#include "StellarEvolution/MESAIO.h"

extern "C" {
    StellarEvolution::MESA::Interpolator* create_interpolator()
    {
        return new StellarEvolution::MESA::Interpolator(
            "../poet_src/StellarEvolution/MESA"
        );
    }
    const StellarEvolution::EvolvingStellarQuantity* 
        interpolate_moment_of_inertia(
            StellarEvolution::MESA::Interpolator* interpolator,
            double mass,
            double metallicity,
            int zone
        )
        {
            return interpolator->interpolate_moment_of_inertia(
                mass,
                metallicity,
                static_cast<Core::StellarZone>(zone)
            );
        }

    double evaluate_quantity(
        StellarEvolution::EvolvingStellarQuantity* quantity,
        double age
    )
    {
        quantity->select_interpolation_region(age);
        return (*quantity)(age);
    }
}
