/**\file
 *
 * \brief The definitions of the functions declared in CInterface.h.
 *
 * \ingroup StellarEvolution
 */

#include "CInterface.h"

const int RADIUS = StellarEvolution::RADIUS;
const int ICONV = StellarEvolution::ICONV;
const int LUM = StellarEvolution::LUM;
const int IRAD = StellarEvolution::IRAD;
const int MRAD = StellarEvolution::MRAD;
const int RRAD = StellarEvolution::RRAD;
const int NUM_QUANTITIES = StellarEvolution::NUM_QUANTITIES;

MESAInterpolator* create_interpolator(const char *mesa_dir,
                                      double *smoothing,
                                      int *nodes,
                                      bool *vs_log_age,
                                      bool *log_quantity)
{
    StellarEvolution::MESA::Interpolator *result;
    if(!smoothing) {
        assert(!nodes);
        result = new StellarEvolution::MESA::Interpolator(mesa_dir);
    } else {
        assert(nodes);
        result = new StellarEvolution::MESA::Interpolator(
            mesa_dir,
            std::vector<double>(smoothing, smoothing + NUM_QUANTITIES),
            std::vector<int>(nodes, nodes + NUM_QUANTITIES),
            std::vector<bool>(vs_log_age, vs_log_age + NUM_QUANTITIES),
            std::vector<bool>(log_quantity, log_quantity + NUM_QUANTITIES)
        );
    }

    return reinterpret_cast<MESAInterpolator*>(result);
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

void evaluate_quantity_array(const EvolvingStellarQuantity *quantity,
                             double *age,
                             unsigned nvalues,
                             double *result)
{
    const StellarEvolution::EvolvingStellarQuantity* actual_quantity =
        reinterpret_cast<const StellarEvolution::EvolvingStellarQuantity*>(
            quantity
        );

    for(unsigned i = 0; i < nvalues; ++i) {
        if(
            i==0
            ||
            age[i] < age[i - 1]
            ||
            actual_quantity->next_discontinuity() < age[i]
        )
            actual_quantity->select_interpolation_region(age[i]);
        result[i] = (*actual_quantity)(age[i]);
    }
}

void differentiate_quantity(const EvolvingStellarQuantity* quantity,
                            double age, double *result)
{
    const StellarEvolution::EvolvingStellarQuantity* actual_quantity =
        reinterpret_cast<const StellarEvolution::EvolvingStellarQuantity*>(
            quantity
        );
    actual_quantity->select_interpolation_region(age);
    const Core::FunctionDerivatives *deriv = actual_quantity->deriv(age);

    for(unsigned order = 0; order < 3; ++order)
        result[order] = deriv->order(order);

    delete deriv;
}

void differentiate_quantity_array(const EvolvingStellarQuantity *quantity,
                                  double *age,
                                  unsigned nvalues,
                                  double *result)
{
    const StellarEvolution::EvolvingStellarQuantity* actual_quantity =
        reinterpret_cast<const StellarEvolution::EvolvingStellarQuantity*>(
            quantity
        );

    for(unsigned i = 0; i < nvalues; ++i) {
        if(
            i==0
            ||
            age[i] < age[i - 1]
            ||
            actual_quantity->next_discontinuity() < age[i]
        )
            actual_quantity->select_interpolation_region(age[i]);
        const Core::FunctionDerivatives 
            *deriv = actual_quantity->deriv(age[i]);
        for(unsigned order = 0; order < 3; ++order)
            result[order * nvalues + i] = deriv->order(order);
        delete deriv;
    }
}

double quantity_min_age(const EvolvingStellarQuantity* quantity)
{
    return 
        reinterpret_cast<const StellarEvolution::EvolvingStellarQuantity*>(
            quantity
        )->range_low();
}

double quantity_max_age(const EvolvingStellarQuantity* quantity)
{
    return 
        reinterpret_cast<const StellarEvolution::EvolvingStellarQuantity*>(
            quantity
        )->range_high();
}

void quantity_continuous_range(const EvolvingStellarQuantity* quantity,
                               double age,
                               double *range_min,
                               double *range_max)
{
    const StellarEvolution::EvolvingStellarQuantity* actual_quantity =
        reinterpret_cast<const StellarEvolution::EvolvingStellarQuantity*>(
            quantity
        );
    actual_quantity->select_interpolation_region(age);
    *range_min = actual_quantity->previous_discontinuity();
    *range_max = actual_quantity->next_discontinuity();
}

void save_interpolator(MESAInterpolator *interpolator, const char *filename)
{
    reinterpret_cast<const StellarEvolution::MESA::Interpolator*>(
        interpolator
    )->save_state(filename);
}

MESAInterpolator *load_interpolator(const char *filename)
{
    StellarEvolution::MESA::Interpolator *interpolator = 
        new StellarEvolution::MESA::Interpolator();
    interpolator->load_state(filename);
    return reinterpret_cast<MESAInterpolator*>(interpolator);
}

double default_smoothing(int quantityID)
{
    assert(quantityID >= 0 && quantityID < NUM_QUANTITIES);
    return StellarEvolution::MESA::Interpolator::default_smoothing(
        static_cast<StellarEvolution::QuantityID>(quantityID)
    );
}

int default_nodes(int quantityID)
{
    assert(quantityID >= 0 && quantityID < NUM_QUANTITIES);
    return StellarEvolution::MESA::Interpolator::default_nodes(
        static_cast<StellarEvolution::QuantityID>(quantityID)
    );
}

bool default_vs_log_age(int quantityID)
{
    assert(quantityID >= 0 && quantityID < NUM_QUANTITIES);
    return StellarEvolution::MESA::Interpolator::default_vs_log_age(
        static_cast<StellarEvolution::QuantityID>(quantityID)
    );
}

bool default_log_quantity(int quantityID)
{
    assert(quantityID >= 0 && quantityID < NUM_QUANTITIES);
    return StellarEvolution::MESA::Interpolator::default_log_quantity(
        static_cast<StellarEvolution::QuantityID>(quantityID)
    );
}
