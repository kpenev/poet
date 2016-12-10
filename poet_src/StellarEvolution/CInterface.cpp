#include "CInterface.h"

const int RADIUS = StellarEvolution::RADIUS;
const int ICONV = StellarEvolution::ICONV;
const int LUM = StellarEvolution::LUM;
const int IRAD = StellarEvolution::IRAD;
const int MRAD = StellarEvolution::MRAD;
const int RRAD = StellarEvolution::RRAD;
const int NUM_QUANTITIES = StellarEvolution::NUM_QUANTITIES;

MESAInterpolator* create_interpolator(const char *mesa_dir)
{
    return reinterpret_cast<MESAInterpolator*>(
        new StellarEvolution::MESA::Interpolator(mesa_dir)
    );
}

void destroy_interpolator(MESAInterpolator *interpolator)
{
    delete reinterpret_cast<StellarEvolution::MESA::Interpolator*>(
        interpolator
    );
}

const EvolvingStellarQuantity* create_quantity(
    const MESAInterpolator* interpolator,
    int quantityID,
    double mass,
    double metallicity
)
{
    assert(quantityID >= 0 && quantityID < NUM_QUANTITIES);
    const StellarEvolution::MESA::Interpolator *actual_interpolator =
        reinterpret_cast<const StellarEvolution::MESA::Interpolator*>(
            interpolator
        );
    return 
        reinterpret_cast<EvolvingStellarQuantity*> (
            (*actual_interpolator)(
                static_cast<StellarEvolution::QuantityID>(quantityID),
                mass,
                metallicity
            )
        );
}

void destroy_quantity(EvolvingStellarQuantity *quantity)
{
    delete reinterpret_cast<StellarEvolution::EvolvingStellarQuantity*>(
        quantity
    );
}

double evaluate_quantity(const EvolvingStellarQuantity* quantity,
                         double age)
{
    const StellarEvolution::EvolvingStellarQuantity* actual_quantity =
        reinterpret_cast<const StellarEvolution::EvolvingStellarQuantity*>(
            quantity
        );
    actual_quantity->select_interpolation_region(age);
    return (*actual_quantity)(age);
}
