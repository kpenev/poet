#include "MESAIO.h"

extern "C" {
    ///Identifier for the stellar radius as an interpolation quantity.
    extern const int RADIUS;

    ///\brief Identifier for the convective zone moment of inertia as an
    ///interpolation quantity.
    extern const int ICONV;

    ///Identifier for the stellar luminosity as an interpolation quantity.
    extern const int LUM;

    ///\brief Identifier for the radiative zone moment of inertia as an
    ///interpolation quantity.
    extern const int IRAD;

    ///\brief Identifier for the radiative zone mass of inertia as an
    ///interpolation quantity.
    extern const int MRAD;

    ///\brief Identifier for the convective-radiative boundary as an
    ///interpolation quantity.
    extern const int RRAD;

    ///The number of interpolation quantities currentyl supported.
    extern const int NUM_QUANTITIES;

    ///Opaque struct to cast to/from StellarEvolution::Interpolator pointers.
    struct MESAInterpolator;

    ///\brief Opaque struct to cast to/from
    ///StellarEvolution::EvolvingStellarQuantity pointers.
    struct EvolvingStellarQuantity;

    ///\brief Create an interpolator from a directory containing MESA tracks.
    ///
    ///The result must be de-allocated when no longer necessary.
    MESAInterpolator* create_interpolator(const char *mesa_dir);

    ///\brief Destroy a previously created interpolator.
    void destroy_interpolator(
        ///The interpolator to destroy. Must have previously been created
        ///using create_interpolator()
        MESAInterpolator *interpolator
    );

    ///\brief Create a single quantity interpolation for a given star.
    ///
    ///The result must be de-allocated when no longer needed.
    const EvolvingStellarQuantity* create_quantity(
        ///An interpolator previously created with create_interpolator().
        const MESAInterpolator* interpolator,

        ///The quantity to interpolate (one of RADIUS, ICONV, LUM, IRAD,
        ///MRAD or RRAD)
        int quantityID,

        ///The stellar mass for which to interpolate the quantity.
        double mass,

        ///The stellar metallicity for which to interpolate the quantity.
        double metallicity
    );

    ///Destroy a previously created evolving stellar quantity.
    void destroy_quantity(
        ///The quantity to destroy. Must have previously been created using
        ///create_quantity().
        EvolvingStellarQuantity *quantity
    );

    ///Evaluate a stellar quantity at a given age.
    double evaluate_quantity(
        ///The quantity to evaluate. Must be previously created using
        ///interpolate().
        const EvolvingStellarQuantity* quantity,

        ///The age at which to evaluate the quantity in Gyrs.
        double age
    );
}
