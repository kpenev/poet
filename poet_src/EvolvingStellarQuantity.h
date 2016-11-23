/**\file
 *
 * \brief Defines classes and methods implementing the intepolation of a
 * single stellar quantity from stellar evolution tracks.
 * 
 * \ingroup StellarSystem_group
 */

#ifndef __EVOLVING_STELLAR_QUANTITY_H
#define __EVOLVING_STELLAR_QUANTITY_H

#include "Functions.h"
#include "Error.h"
#include <valarray>
#include <list>
#include <string>
#include <iostream>
#include <fstream>

///Perform a bi-cubic spline interpolation of a single quantity.
double mass_metallicity_interp(
    ///The masses of the stelar models on which to base the interpolation in 
    ///\f$M_\odot\f$
    alglib::real_1d_array interp_masses,

    ///The metallicities of the stellar models on which to base the
    ///interpolation in \f$[Fe/H]\f$.
    alglib::real_1d_array interp_metallicities,

    ///The values of the quantity being interpolated on the grid defined by
    ///\p interp_masses and \p interp_metallicities.
    alglib::real_1d_array interp_values,

    ///The stellar mass to which to interpolate in \f$M_\odot\f$.
    double stellar_mass,

    ///The stellar metallicity to which to interpolate in \f$[Fe/H]\f$.
    double stellar_metallicity
);

///\brief A class that calculates derivatives with respect to an argument for
///functions of the log(argument).
///
///\ingroup StellarSystem_group
class LogArgDerivatives : public FunctionDerivatives {
private:
	///The value of the argument at which derivatives are calculated.
	double x;

	///The value (zeroth derivative).
	mutable double value;

	///\brief Was the interpolation done against the logarithm of the
	///argument and hence the derivative needs to be corrcted.
	bool correct_log_arg;

	///\brief The currently computed derivatives.
	///
	///With respect to log(arg) if correct_log_arg is true, or with respect
	///to arg if not.
	mutable std::vector<double> underlying_deriv_values,

			///Previously calculated values of the derivatives.
			///
			///These are reused if requested multiple times.
			deriv_values;

	///\brief Actually corrects for differentiating w.r.t. log(arg) instead
	///of arg.
	///
	///Uses the pre-computed array of the derivatives up to the given order
	///in underlying_deriv_values with respect to ln(x) to return the
	///order-th derivative with respect to x.
	double transform_log_arg_deriv(unsigned order) const;
protected:
	///\brief Should be overwritten to calculate the derivatives with respect
	///to either arg or log(arg) as specified on construction.
	virtual double calc_deriv(unsigned deriv_order) const =0;
public:
	///\brief Create a derivative for functions of possibly log(arg).
	///
	///The created object corrects for the fact that the underlying
	///derivative (defined by the calc_deriv method) is with respect to the
	///logarithm of the argument if arg_val is not NaN. No correction if it
	///is Nan.
	LogArgDerivatives(double arg_val=NaN) :
		x(arg_val), value(NaN), correct_log_arg(!std::isnan(arg_val)) {}

	///Returns the deriv_order-th derivative of the quantity
	double order(unsigned deriv_order=1) const;
};

///\brief Derivative class for stellar quantities which are interpolated
///age, mass and metallicity.
///
///\ingroup StellarSystem_group
class InterpolatedDerivatives : public LogArgDerivatives {
private:
	double 
        ///The mass to interpolate to in \f$M_\odot\f$.
        __stellar_mass,
        
        ///The metallicity to interpolate to in Solar metallicities;
        __stellar_metallicity;

	///The age derivatives for each stellar model.
	std::vector<const FunctionDerivatives *> *__interp_deriv;

	
    const alglib::real_1d_array 
        ///The masses of the stelar models in \f$M_\odot\f$
        &__interp_masses,

        ///The metallicities of the stellar models in Solar metallicities.
        &__interp_metallicities;
        
	///Whether to delete the derivatives it was created with
	bool __delete_derivatives;
protected:
	///Returns the deriv_order-th derivative of the quantity
	double calc_deriv(unsigned deriv_order) const;
public:
	///\brief Create an object that interpolates derivatives from evolution
	///tracks.
	///
	///If age is specified the input derivatives are assumed to be with
	///respect to ln(age), while derivatives always with respect to age are
	///output.
	InterpolatedDerivatives(
        double mass,
        double metallicity,
        std::vector<const FunctionDerivatives*> *derivatives,
        const alglib::real_1d_array &interp_masses,
        const alglib::real_1d_array &interp_metallicities,
        double age = NaN,
        bool delete_derivatives = false
    );
	
	///Deletes the interpolation data if so specified on creation.
	~InterpolatedDerivatives()
	{
		if(delete_derivatives) {
			for(size_t i=0; i<interp_deriv->size(); i++) 
				delete (*__interp_deriv)[i];
			delete __interp_deriv; 
		}
	}
};

///\brief Makes a derivative with respect to linear argument from a
///derivative with respect to log(argument).
class RemoveLogDeriv : public LogArgDerivatives {
private:
	///The original logarithmic derivative 
	const FunctionDerivatives *__log_deriv;

	///Whether to delete the underlying log-derivative on destruction.
	bool __delete_deriv;
protected:
	///Returns the deriv_order-th derivative of the quantity
	double calc_deriv(unsigned deriv_order) const
	{return __log_deriv->order(deriv_order);}
public:
	///Create a linear derivative from a log one.
	RemoveLogDeriv(double age, const FunctionDerivatives *log_deriv,
			bool delete_deriv) :
		LogArgDerivatives(age), __log_deriv(log_deriv), 
		__delete_deriv(delete_deriv) {}

	///Deletes the input logarithmic derivative if so specified on creation.
	~RemoveLogDeriv()
	{if(__delete_deriv) delete __log_deriv;}
};

///\brief Derivative class for stellar quantities which age scaled
///quantities for a tabulated mass.
///
///\ingroup StellarSystem_group
class ScaledDerivatives : public LogArgDerivatives {
private:
	///The derivative of the underlying quantity
	const FunctionDerivatives *underlying_deriv;

	///The scaling applied to the age of the underlying quantity
	double scaling;

	///Whether to delete the input derivative when the object is destroyed.
	bool delete_underlying,
		 
		 ///\brief Was the interpolation done against log(argument) and hence
		 ///the derivative needs to be corrcted.
		 correct_log_arg;
protected:
	///\brief Returns the deriv_order-th derivative of the quantity.
	///
	///It should return the derivative with respect to either log(arg) or arg
	///as specified on construction.
	double calc_deriv(unsigned deriv_order) const;
public:
	///\brief Construct a derivative from the derivative of a function whose
	///argument is scaled by the given factor.
	///
	///If age is given and not NaN,
	///the underlying derivatives are assumed to be with respect to ln(age)
	///and a correction is made to return derivatives with respect to age.
	ScaledDerivatives(const FunctionDerivatives *deriv, double factor,
			double age=NaN, bool delete_deriv=false) :
		LogArgDerivatives(age), underlying_deriv(deriv), scaling(factor),
		delete_underlying(delete_deriv), correct_log_arg(!std::isnan(age)) {}

	///Deletes the derivative if so specified on creation.
	~ScaledDerivatives() {
		if(delete_underlying) delete underlying_deriv;
	}

};

///\brief Derivative class for a quantity that is the sum of two other
///quantities.
///
///\ingroup StellarSystem_group
class SumDerivatives : public FunctionDerivatives {
private:
	///The derivatives of the first quantity in the sum.
	const FunctionDerivatives *q1_deriv,

		  ///The derivatives of the second quantity in the sum.  
		  *q2_deriv;

	///Whether to delete the input derivative when the object is destroyed.
	bool destroy_derivs;
public:
	///Create a derivative object for a sum of two quantities: q1+q2.
	SumDerivatives(
			///Pointer to the derivative of the first quantity (q1).
			const FunctionDerivatives *derivative1,

			///Pointer to the derivative of the second quantity (q2).
			const FunctionDerivatives *derivative2,
			
			///Delete the input derivatives on destruction?
			bool delete_inputs=true)
		: q1_deriv(derivative1), q2_deriv(derivative2),
		destroy_derivs(delete_inputs) {}
	
	///The deriv_order-th derivative.
	double order(unsigned deriv_order=1) const
	{return q1_deriv->order(deriv_order)+q2_deriv->order(deriv_order);}

	///Clean up.
	~SumDerivatives()
	{if(destroy_derivs) {delete q1_deriv; delete q2_deriv;}}
};

///The four directions in which an interpolation grid can grow.
class AllowedGridGrowth {
private:
    bool 
        ///Can grow toward lower masses?
        __lighter,

        ///Can grow toward higher masses?
        __heavier,

        ///Can grow toward lower metallicities?
        __poorer,

        ///Can grow toward higher metallicities?
        __richer;

public:
    AllowedGridGrowth()
        : __lighter(true), __heavier(true), __poorer(true), __richer(true)
    {}

    ///Can grow toward lower masses?
    bool lighter() const {return __lighter;}

    ///Can grow toward higher masses?
    bool heavier() const {return __heavier;}

    ///Can grow toward lower metallicities?
    bool poorer() const {return __poorer;}

    ///Can grow toward higher metallicities?
    bool richer() const {return __richer;}

    ///Disable growth to lower masses.
    AllowedGridGrowth &block_lighter() {__lighter = false; return *this;}

    ///Disable growth to higher masses.
    AllowedGridGrowth &block_heavier() {__heavier = false; return *this;}

    ///Disable growth to lower metallicities.
    AllowedGridGrowth &block_poorer() {__poorer = false; return *this;}

    ///Disable growth to higher metallicities.
    AllowedGridGrowth &block_richer() {__richer = false; return *this;}

    ///Is growth allowed in at least one direction?
    operator bool() const
    {return __lighter || __heavier || __poorer || __richer;}
}

///\brief A class for stellar properties that depend on age.
///
///\ingroup StellarSystem_group
class EvolvingStellarQuantity : public OneArgumentDiffFunction {
private:
    double 
        ///The mass to which to interpolate in \f$M_\odot\f$.
        __mass,

        ///The mass to which to interpolate in Solar metallicities.
        __metallicity;

	
	double 
        ///The minimum age for which this quantity is defined in Gyr.
        __min_age,

        ///The maximum age for which this quantity is defined in Gyr.
        __max_age; 

	///Whether the tracks have log(age) instead of age as their argument.
	bool __use_log_age,
		 
		 ///Should the quantity be assumed zero below the minimum track age.
		 __initially_zero;
	
	const std::valarray<double> 
        ///The masses of the evolution tracks below the high low mass split.
        &__track_masses,

        ///The metallicities of the evolution tracks.
        &__track_metallicities;

    std::valarray<double> 
        ///\brief The minimum interpolation age for the current star to which 
        ///each track can contribute.
        __min_interp_ages,

        ///\brief The maximum interpolation age for the current star to which 
        ///each track can contribute.
        __max_interp_ages;

    ///The ages at which the interpolation grid needs to be re-determined.
    std::vector<double> __interp_grid_change_ages;

    ///\brief The entry in ::__interp_grid_change_ages up to which the
    ///current interpolation grid is valid.
    std::vector<double>::const_iterator __next_grid_change_age;
    
    size_t 
        ///The index of the smallest track mass not smaller than ::__mass.
        __mass_index_above,

        ///The index of the largest track mass not exceeding ::__mass.
        __mass_index_below,

        ///\brief The index of the smallest track metallicity not smaller
        ///than ::__metallicity.
        __metallicity_index_above,

        ///\brief The index of the largest track metallicity not exceeding
        ///::__metallicity.
        __metallicity_index_below;

	///\brief The model tracks for the evolution of the quantity on the grid
    ///defined by ::__track_masses and ::__track_metallicities.
    ///
    ///The mass index varies faster
	std::vector<const OneArgumentDiffFunction *> __evolution_tracks;

    ///\brief How far to extrapolate past the last tabulated age for
    ///a model.
    ///
    ///Models are included in the mass interpolation only if
    ///the required age is no larger than the maximum tabulated age 
    ///times this factor.
    double __extrapolate;

    alglib::real_1d_array 
        ///The current track masses participating in the interpolation.
        __interp_masses,

        ///\brief The current track metallicities participating in the
        ///interpolation.
        __interp_metallicities;

    size_t 
        ///\brief The index within ::__track_masses of the lowest mass
        ///currently participating in the interpolation.
        __min_interp_mass_index,

        ///\brief The index within ::__track_masses of the highest mass
        ///currently participating in the interpolation.
        __max_interp_mass_index,

        ///\brief The index within ::__track_metallicities of the lowest
        ///metallicity currently participating in the interpolation.
        __min_interp_metallicity_index,

        ///\brief The index within ::__track_metallicities of the highest
        ///metallicity currently participating in the interpolation.
        __max_interp_metallicity_index;

    ///\brief Return the index within ::__evolution_tracks for the given mass
    ///and metallicity indices.
    inline const size_t &track_index(
        ///The index within ::__track_masses of the desired mass.
        size_t mass_index,

        ///The index within ::__track_metallicities of the desired
        ///metallicity.
        size_t metallicity_index) const
    {return metallicity_index * __track_masses.size() + mass_index;}

    ///\brief Answer if a given track can participate in interpolating to the
    ///given age.
    inline bool track_in_range(
        ///The index of the track to check within ::__evolution_tracks.
        size_t track_i,

        ///The age to which interpolation is desired.
        double age) const
    {return (__min_interp_ages[track_i] <= age
             &&
             __max_interp_ages[track_i] >= age);}

    ///\brief Answer if a given track can participate in interpolating to the
    ///given age.
    inline bool track_in_range(
        ///The index of the mass of the track to check within
        ///::__track_masses.
        size_t mass_i,

        ///The index of the metallicity of the track to check within
        ///::__track_metallicities.
        size_t metallicity_i,

        ///The age to which interpolation is desired.
        double age) const
    {return track_in_range(track_index(mass_i, metallicity_i));}

    ///\brief Verify that the stellar mass and metallicity are within range
    ///of the evolution tracks.
    void check_grid_range() const;

    ///\brief The two indices within the given sorted array defining the
    ///closed internal containing value.
    ///
    ///If value is not exactly equal to an array entry, the two indices are
    ///consecutive, if the value is exactly equal to an entry, the two
    ///indices are the same.
    void find_cell(
        ///The boundaries of the grid cells in a single dimension.
        const std::valarray<double> &boundaries

        ///The value whose cell we are looking for.
        const double &value,

        ///The index of the largest boundary <= \p value.
        size_t &below_index,

        ///The index of the smallest boundary >= \p value.
        size_t &above_index,
    );

    ///Fill the ::__min_interp_ages and ::__max_interp_ages members.
    void set_interp_age_ranges();

	///\brief Interpolate the quantity for the given track to the given age,
	///returning NaN if out of age range.
	///
	///If derivatives is not NULL initializes that to a pointer to a
	///derivatives at the current age structure.
	double evaluate_track(
        ///The age at which to evaluate the track. If ::__use_log_age is
        ///true, the track is evaluated at log(\p track_age), otherwise it is
        ///directly passed to the track. In particular, the caller should
        ///have already transformed this argument to the correct
        ///interpolation parameter.
        double track_age,

        ///The track to evaluate. Should be an entry from
        ///::__evolution_tracks.
        const OneArgumentDiffFunction &track,

        ///If not NULL, \p *derivatives is set to a newly allocated
        ///derivative instance of one of the children FunctionDerivatives
        ///classes.
        const FunctionDerivatives **derivatives) const;

    ///\brief Figure out in which directions we can expand a mass-metallicity
    ///interpolation grid, assuming a single direction expansion.
    void check_grid_expansion_directions(
        ///The current state of the grid expansion possibilities. On output,
        ///directions in which growth is no longer allowed are disabled.
        AllowedGridGrowth &grow
    ) const;

    ///Expand the current interpolation grid in one of the allowed
    ///directions.
    void expand_grid(
        ///The current state of the grid expansion possibilities.
        const AllowedGridGrowth &grow
    ) const;

    ///Find the best sub-grid of tracks to interpolate on.
    void update_interpolation_range(
        ///The age to which we are trying to interpolate.
        double age
    );

	///\brief Interpolate the quantity to the desired age assuming we are in
	///the low mass regime.
	///
	///If derivatives is not NULL initializes that to a pointer to a
	///derivatives at the current age structure.
	double interpolate(double age, 
                       const FunctionDerivatives **derivatives=NULL) const;

    ///\brief Return the interpoltaion parameter for the given age for the
    ///current star.
    template<typename VALUE_TYPE>
    inline const VALUE_TYPE &age_to_interp_param(
        const VALUE_TYPE &age
    )
    {return age_to_interp_param(age, __mass, __metallicity);}

    ///\brief Return the age for the given interpoltaion parameter for the
    ///current star.
    template<typename VALUE_TYPE>
    inline const VALUE_TYPE &interp_param_to_age(
        const VALUE_TYPE &interp_param
    )
    {return interp_param_to_age(interp_param, __mass, __metallicity);}

protected:
    ///\brief Return the interpoltaion parameter given age, mass and
    ///metallicity.
    ///
    ///Must be an increasing monotonic function
    template<typename VALUE_TYPE>
    virtual const VALUE_TYPE &age_to_interp_param(
        ///The age for which the interpolation parameter is needed in Gyrs.
        const VALUE_TYPE &age,
        
        ///The stellar mass for which the interpolation parameter is needed
        ///in \f$M_\odot\f$.
        const VALUE_TYPE &mass,

        ///The stellar metallicity for which the interpolation parameter is
        ///needed in \f$[Fe/H]\f$.
        const VALUE_TYPE &metallicity);

    ///\brief Return the age in Gyrs given an interpolation parameter, mass,
    ///and metallicity.
    ///
    ///Must be an increasing monotonic function
    template<typename VALUE_TYPE>
    virtual const VALUE_TYPE &interp_param_to_age(
        ///The interpolation parameter for which the age is needed.
        const VALUE_TYPE &interp_param,
        
        ///The stellar mass for which the interpolation parameter is needed
        ///in \f$M_\odot\f$.
        const VALUE_TYPE &mass,

        ///The stellar metallicity for which the interpolation parameter is
        ///needed in \f$[Fe/H]\f$.
        const VALUE_TYPE &metallicity);

public:
	///\brief Construct an object that can be set to interpolate between
	///tabulated evolution tracks of a quantity.
	EvolvingStellarQuantity() {};

	///Create an evolving quantity that interpolates to the given mass.
	EvolvingStellarQuantity(
			///The stellar mass to interpolate to in \f$M_\odot\f$
			double mass, 

			///The stellar metallicity (\f$[Fe/H]f$) to interpolate to.
			double metallicity, 

			///The masses for which evolution tracks are given in
			/// \f$M_\odot\f$
			const std::valarray<double> &track_masses,

			///The metallicities (\f$[Fe/H]\f$) for which evolution tracks
            ///are given.
			const std::valarray<double> &track_metallicities,

			///The evolution tracks of the relevant quantity on the grid
            ///defined by \p track_masses and \p track_metallicities. The
            //mass index varies faster.
			const std::list<const OneArgumentDiffFunction *> 
            &evolution_tracks,

			///Whether the track uses log(age) as the independent argument
			///instead of age.
			bool log_age=true,

			///Model tracks are included in the mass-metallicity
            ///interpolation only if the required age is no larger than the
            ///maximum tabulated age times this factor.
			double extrapolate=1.01,

			///Whether this is a quantity that is identically zero below some
			///age and turns on afterwards
			bool starts_zero=false);

    ///\brief Prepare the quantity for interpolation around the given age.
    ///
    ///After calling this method, requesting values or derivatives outside
    ///the range of the continuous region containing this age
    ///(see ::discontinuities) fails an assert.
    virtual void select_interpolation_region(double age);

	///Return the value the quantity takes at the given age in Gyr.
	virtual double operator()(double age) const
    {return interpolate(age);}

	///Return the age derivative of the quantity at the given age in Gyr.
	virtual const FunctionDerivatives *deriv(double age) const;

	///The largest age for which the quantity can be interpolated in Gyr.
	virtual double range_high() const {return max_age;}

	///The smallest age for which the quantity can be interpolated in Gyr.
	virtual double range_low() const {return min_age;}

    ///The ages at which the quantity may be discontinuous.
    virtual const std::vector<double> &discontinuities() const
    {return __interp_grid_change_ages;}

    ///\brief The upper bound of the current interpolation region (over which 
    ///the quantity is guaranteed continuous).
    ///
    ///The current interpolation region is set either
    ///by ::set_starting_interpolation_age or
    ///by ::enable_next_interpolation_region.
    virtual double next_discontinuity() const
    {return *__next_grid_change_age;}

    ///\brief Set up the interpolation over the next interpolation region
    ///(between consecutive discontinuities.)
    virtual double enable_next_interpolation_region();

	///An iterator over the ages (in Gyr) where the quantity takes the given
	///y value.
	InterpSolutionIterator crossings(double =0) const
	{throw Error::Runtime("Called EvolvingStellarQuantity::crossings, "
			"which are ill defined.");}
};

class ZeroQuantity : public EvolvingStellarQuantity {
public:
    ///Do nothing. See EvolvingStellarQuantity::select_interpolation_region.
    void select_interpolation_region(double age) {}

	///Return the value the quantity takes at the given age.
	double operator()(double) const {return 0;}

	///Return the age derivative of the quantity at the given age.
	const FunctionDerivatives *deriv(double) const
	{return new ZeroDerivatives;}

	///The largest age for which the quantity can be interpolated
	double range_high() const {return Inf;}

	///The smallest age for which the quantity can be interpolated.
	double range_low() const {return -Inf;}

    ///No discontinuities. See EvolvingStellarQuantity::next_discontinuity.
    double next_discontinuity() const {return Inf;}

    ///\brief Do nothing.
    ///See EvolvingStellarQuantity::enable_next_interpolation_region.
    double enable_next_interpolation_region();

	///An iterator over the ages where the quantity takes the given y value.
	InterpSolutionIterator crossings(double =0) const
	{throw Error::Runtime("Called ZeroQuantity::crossings, "
			"which are ill defined.");}
};

///\brief A clas for stellar quantities that are the sum of two other
///quantities.
///
///\ingroup StellarSystem_group
class SumQuantity : public EvolvingStellarQuantity {
private:
	///This quantity will be q1+q2
	const EvolvingStellarQuantity *q1, *q2;

	///Whether to destroy the input quantities on destruction
	bool destroy_qs;
public:
	///Create a quantity that is (*quantity1)-(*quantity2)
	SumQuantity(const EvolvingStellarQuantity *quantity1,
			const EvolvingStellarQuantity *quantity2,
			bool delete_inputs=false)
		: q1(quantity1), q2(quantity2), destroy_qs(delete_inputs) {}

    ///See EvolvingStellarQuantity::select_interpolation_region.
    virtual void select_interpolation_region(double age)
    {
        q1->select_interpolation_region(age);
        q2->select_interpolation_region(age);
    }

	///Return the value the quantity takes at the given age.
	double operator()(double age) const
	{return (*q1)(age)+(*q2)(age);}

	///Return the age derivative of the quantity at the given age.
	const FunctionDerivatives *deriv(double age) const
	{return new SumDerivatives(q1->deriv(age), q2->deriv(age), true);}

	///The largest age for which the quantity can be interpolated
	double range_high() const
	{return std::min(q1->range_high(), q2->range_high());}

	///The smallest age for which the quantity can be interpolated.
	double range_low() const
	{return std::max(q1->range_low(), q2->range_low());}

    ///See EvolvingStellarQuantity::next_discontinuity.
    double next_discontinuity() const
    {return std::min(q1->next_discontinuity(), q2->next_discontinuity());}

    ///See EvolvingStellarQuantity::enable_next_interpolation_region.
    double enable_next_interpolation_region()
    {
        if(q1->next_discontinuity() == q2->next_discontinuity()) {
            q1->enable_next_interpolation_region();
            q2->enable_next_interpolation_region();
        } else if(q1->next_discontinuity() < q2->next_discontinuity())
            q1->enable_next_interpolation_region();
        else 
            q2->enable_next_interpolation_region();
    }

	///An iterator over the ages where the quantity takes the given y value.
	InterpSolutionIterator crossings(double =0) const
	{throw Error::Runtime("Called EvolvingStellarQuantity::crossings, "
			"which are ill defined.");}

	///Clean up.
	~SumQuantity()
	{if(destroy_qs) {delete q1; delete q2;}}
};

#endif
