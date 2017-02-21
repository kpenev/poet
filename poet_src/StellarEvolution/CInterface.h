/**\file
 *
 * \brief Declare C-style functions for accessing the functionality of the
 * StellarEvolution library.
 *
 * \ingroup StellarEvolution
 */

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
    MESAInterpolator* create_interpolator(
        ///The directory containing all and only the MESA tracks to include
        ///in the interpolation.
        const char *mesa_dir,
        
        ///The set of smoothing arguments to use. One for each quantity, in
        ///the order defined by the RADIUS, ICONV, ..., constants.
        double *smoothing,

        ///The set of interpolation nodes to use in the same order as
        ///smoothing.
        int *nodes,

        ///For each quantity the corresponding entry decides if log(age) will
        ///be the independent argument against which interpolation is
        ///performed (instead of age).
        bool *vs_log_age,

        ///For each quantity the corresponding entry decides if log(quantity)
        ///will be interpolated instead of just quantity.
        bool *log_quantity
    );

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

    ///Evaluate a stellar quantity at an array of ages.
    void evaluate_quantity_array(
        ///The quantity to evaluate. Must be previously created using
        ///interpolate()
        const EvolvingStellarQuantity *quantity,
        
        ///The array of ages to evaluate the quantity at.
        double *age,

        ///The number of ages at which evaluation is required.
        unsigned nvalues,

        ///A pre-allocated memory (size: nvalues) where to place the result.
        double *result
    );

    ///Calculate the zeroth, first and second derivatives of a quantity.
    void differentiate_quantity(
        ///The quantity to differentiate. Must be previously created using
        ///interpolate().
        const EvolvingStellarQuantity* quantity,

        ///The age at which to differentiate the quantity in Gyrs.
        double age,

        ///A pre-allocated array of size 3 where to place the result.
        double *result
    );

    ///\brief Calculate the derivatives of a quantity at an array of ages.
    ///
    ///The result is a double array containing 3 sub-arrays of
    ///zeroth, first and second order derivatives. That is the first
    ///consecutive nvalues entries are the function values at echo of the
    ///nvalues ages, the next nvalues entries are the first derivatives at
    ///each age etc. 
    void differentiate_quantity_array(
        ///The quantity to differentiate. Must be previously created using
        ///interpolate()
        const EvolvingStellarQuantity *quantity,

        ///The array of ages to evaluate the quantity at.
        double *age,

        ///The number of ages at which evaluation is required.
        unsigned nvalues,

        ///A pre-allocated memory (size: 3 * nvalues) where to place the 
        ///result.
        double *result
    );

    ///Return the minimum age for which the quantity is defined.
    double quantity_min_age(
        ///The quantity to characterize. Must be previously created using
        ///interpolate().
        const EvolvingStellarQuantity* quantity
    );

    ///Return the maximum age for which the quantity is defined.
    double quantity_max_age(
        ///The quantity to characterize. Must be previously created using
        ///interpolate().
        const EvolvingStellarQuantity* quantity
    );

    ///\brief Return the range of ages surrounding a given age over which a
    ///quantity is guaranteed continuous.
    void quantity_continuous_range(const EvolvingStellarQuantity* quantity,
                                   double age,
                                   double *range_min,
                                   double *range_max);

    ///Save the state of an interpolator for faster creation.
    void save_interpolator(
        ///The interpolator to save. Must have previously been created
        ///using create_interpolator()
        MESAInterpolator *interpolator,

        ///The name of the file to save the state to.
        const char *filename
    );

    ///Load a previously saved interpolator state (faster than creating it).
    MESAInterpolator *load_interpolator(
        ///The name of the file to save the state to.
        const char *filename
    );

    ///Return the default smoothing argument used for the given quantity.
    double default_smoothing(
        ///The quantity to return the default smoothing for.
        int quantityID
    );

    ///\brief Return the default number of interpolation nodes used for the
    ///given quantity.
    int default_nodes(
        ///The quantity to return the default nodes for.
        int quantityID
    );

    ///\brief Return whether by default the given quantity is interpolated
    ///vs. log(age).
    bool default_vs_log_age(
        ///The quantity to return the default nodes for.
        int quantityID
    );

    ///\brief Return whether by default the log(given quantity) is interpolated
    ///vs. the quantity itself.
    bool default_log_quantity(
        ///The quantity to return the default nodes for.
        int quantityID
    );
} //End extern "C"
