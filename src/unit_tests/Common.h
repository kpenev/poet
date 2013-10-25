/**\file
 *
 * \brief Functions and classes of general use for all unit tests.
 *
 * \ingroup UnitTests_group
 */

#ifndef __UNITTEST_COMMON_H
#define __UNITTEST_COMMON_H

#include "../Functions.h"
#include "../StellarEvolution.h"
#include "../StellarSystem.h"
#include "../Star.h"
#include "../Planet.h"
#include "../OrbitSolver.h"
#include <cpptest.h>
#include <valarray>
#include <sstream>
#include <gsl/gsl_roots.h>

///\brief Returns true iff \f$|x-y|\leq\mathrm{abs\_tolerance} +
/// \mathrm{frac\_tolerance}\cdot\max(|x|,|y|)\f$.
bool check_diff(double x, double y, double frac_tolerance, 
		double abs_tolerance);

///\brief Returns true iff \f$ \forall i\ |x_i-y_i|\leq
/// \mathrm{abs\_tolerance}_i +
/// \mathrm{frac\_tolerance}_i\cdot\max(|x_i|,|y_i|)\f$.
bool check_diff(std::valarray<double> x, std::valarray<double> y,
		std::valarray<double> frac_tolerance,
		std::valarray<double> abs_tolerance);

///\todo Get rid of this function and use check_diff instead.
bool isEqual(double a, double b);

///\todo Get rid of this function and use check_diff instead.
double getError(double predicted, double actual);

///\todo Get rid of this function and use check_diff instead.
bool approxEqual(double predicted, double actual, double thres=0.02);

///The lowest stellar mass to use in tests in \f$M_\odot\f$.
const double min_stellar_mass=0.4,
	  
	  ///The boundary between high and low mass stars in \f$M_\odot\f$.
	  max_low_mass=1.075,

	  ///The highest stellar mass to use in tests in \f$M_\odot\f$.  
	  max_stellar_mass=1.3,

	  ///Most tests start at this age in Gyr.
	  min_age=1e-7,

	  ///Most tests end at this age in Gyr.  
	  max_age=10.0;

///The lower limit of the mass of random planets.
const double min_planet_mass=10,

	  ///The upper limit of the mass of random planets.
	  max_planet_mass=80,

	  ///The lower limit of the radius of random planets.
	  min_planet_radius=5,

	  ///The upper limit of the radius of random planets.
	  max_planet_radius=15;

///\brief Outputs a comma separated list of the values in the array to the
///given stream.
std::ostream &operator<<(std::ostream &os,
		const std::valarray<double> &array);

///\brief An EvolvingStellar quantity that uses a polynomial instead of
///interpolating.
///
///It also serves as its own derivative.
///
///\ingroup UnitTests_group
class PolynomialEvolutionTrack : public EvolvingStellarQuantity,
	FunctionDerivatives {
private:
	///The coefficients of the polynomial giving the evolution.
	std::valarray<double> poly_coef;

	///The lower limit of the age at which this quantity can be evaluated.
	double xmin,

		   ///\brief The upper limit of the age at which this quantity can be
		   ///evaluated.
		   xmax;

	///The location at which the derivative has been requested.
	mutable	double deriv_x;
public:
	///Create an evolving quantity or its derivative.
	PolynomialEvolutionTrack(

			///The coefficients defining the polynomial that gives the
			///evolution of the quantity.
			const std::valarray<double> &coefficients,

			///The lower end of the range for which the quantity is defined.
			double range_low,

			///The upper end of the range for which the quantity is defined.
			double range_high,
			
			///The abscissa at which the derivatives are to be evaluated.
			double derivative_x=NaN) :
		poly_coef(coefficients), xmin(range_low), xmax(range_high),
		deriv_x(derivative_x) {}

	///\brief Evaluates the polynomial at the given x.
	///
	///Throws an exception if x is outside the range over which the
	///quantity is defined. 
	double operator()(double x) const {deriv_x=x; return order(0);}

	///The lower end of the range where the function is defined
	double range_low() const {return xmin;}

	///The upper end of the range where the function is defined
	double range_high() const {return xmax;}

	///The derivatives.
	const PolynomialEvolutionTrack *deriv(double x) const
	{return new PolynomialEvolutionTrack(poly_coef, xmin, xmax, x);}

	///The order-th derivative
	double order(unsigned deriv_order=1) const;

	///No crossing calculated.
	InterpSolutionIterator crossings(double y=0) const
	{return InterpSolutionIterator();}

	///Prints an expression of the polynomial that this track represents.
	friend std::ostream &operator<<(std::ostream &os,
			const PolynomialEvolutionTrack &track);

	///A string identifying the type of quantity this is.
	std::string kind() const {return "PolynomialEvolutionTrack";}
};

///\brief Implements a StellarEvolution based on polynomial evolution
///quantities.
///
///\ingroup UnitTests_group
class MockStellarEvolution : public StellarEvolution {
private:
	///The polynomial coefficients of the radius dependence on mass and age.
	std::valarray< std::valarray<double> > R,

		///\brief The polynomial coefficients of the convective zone moment
		///of inertia dependence on stellar mass and age.
		Iconv,

		///\brief The polynomial coefficients of the radiative zone moment of
		///inertia dependence on stellar mass and age.
		Irad,

		///\brief The polynomial coefficients of total stellar moment of
		///inertia dependence on stellar mass and age.
		Itot,

		///\brief The polynomial coefficients of the core radius dependence
		///on stellar mass and age.
		Rcore,

		///\brief The polynomial coefficients of the core mass dependence
		///on stellar mass and age.
		Mcore,

		///\brief The polynomial coefficients of the envelope mass dependence
		///on stellar mass and age.
		Menv;

	///The age in Gyr at which the core forms.
	double core_formation;
public:
	///\brief Uses the given values and polynomials for the quantities.
	///
	///If empty polynomial coefficients are specified or NaN for constant
	///quantities in which case random polynomials or values are used.
	MockStellarEvolution(
			///The age at which the core forms in Gyr.
			double core_formation_age=NaN,

			///The polynomial coefficients of the radius dependence on mass
			///and age.
			const std::valarray< std::valarray<double> > &radius=
				std::valarray< std::valarray<double> >(),

			///The polynomial coefficients of the convective zone moment of
			///inertia dependence on stellar mass and age.
			const std::valarray< std::valarray<double> >
				&conv_moment_of_inertia=
				std::valarray< std::valarray<double> >(),

			///The polynomial coefficients of the radiative zone moment of
			///inertia dependence on stellar mass and age.
			const std::valarray< std::valarray<double> > 
				&rad_moment_of_inertia=
				std::valarray< std::valarray<double> >(),

			///The polynomial coefficients of the core radius dependence on
			///stellar mass and age.
			const std::valarray< std::valarray<double> > &core_radius=
				std::valarray< std::valarray<double> >(),

			///The polynomial coefficients of the core mass dependence on
			///stellar mass and age.
			const std::valarray< std::valarray<double> > &core_mass=
				std::valarray< std::valarray<double> >());

	///\brief The Moment of inertia of the specified zone of a star of the
	///specified mass as a function of age.
	///
	///The result must be destroyed when it becomes obsolete.
	const EvolvingStellarQuantity *interpolate_moment_of_inertia(
			double stellar_mass, StellarZone zone, double present_age=-1)
		const;

	///\brief The radius of a  star of the specified mass as a function of
	///age.
	///
	///The result must be destroyed when it becomes obsolete.
	const EvolvingStellarQuantity *interpolate_radius(
			double stellar_mass, double present_age=-1) const;

	///\brief The mass of of the specified zone of a star of the specified
	///mass as a function of age.
	///
	///The result must be destroyed when it becomes obsolete.
	const EvolvingStellarQuantity *interpolate_zone_mass(
			double stellar_mass, StellarZone zone) const;

	///\brief The radius of the core for a star of the specified mass as a
	///function of age.
	///
	///The result must be destroyed when it becomes obsolete.
	const EvolvingStellarQuantity *interpolate_core_boundary(
			double stellar_mass, double present_age=-1) const;

	///The age in Gyr at which the core forms.
	double core_formation_age() const {return core_formation;}
};


///\brief A messy class able to create random stars.
///
///\ingroup UnitTests_group
class StarData {
private:
	unsigned num_stars_created;
public:
	double mass, radius, age, conv_spin, rad_spin, tidal_Q, wind_strength,
		   wind_sat_freq, coupling_timescale, Q_trans_width, disk_lock_w,
		   disk_lock_time;
	PolynomialEvolutionTrack *Lrad_track, *Lconv_track;
	std::valarray<double> evolution_masses, evolution_ages, Lrad, Lconv,
		Iconv, Irad, Mrad_deriv, Lconv_deriv, Lrad_deriv, Rrad, all_radii;
	std::valarray< std::valarray<double> > r_coef, Iconv_coef, Itot_coef,
		Mrad_coef, Rcore_coef;

	StarData();
	
	///Creates a star with random properties (mass, radius, age, ..., random
	///polynomial time functions as moments of inertia ...).
	void create_random_star(Star **star);

	///Deletes the convective and radiative angular momenta tracks if they
	///were allocated.
	~StarData();
};

///\brief A messy class able to create random planets.
///
///\ingroup UnitTests_group
class PlanetData {
public:
	Star* star;
	StarData* sdata;
	double mass, radius;
	std::valarray<double> ages;
	std::valarray<double> semis;
	PlanetData() : star(NULL) {
		sdata = new StarData();
		sdata->create_random_star(&star);
	};
	void create_random_planet(Planet** planet);
	double get_semi(double age);

	///Cleanup.
	~PlanetData();
};

///\brief A messy class able to create random planet-star systems.
///
///\ingroup UnitTests_group
class SystemData {
public:
	Star* star;
	Planet* planet;
	StarData* sdata;
	PlanetData* pdata;
	SystemData() {
		//exit(-1);
		pdata = new PlanetData();
		pdata->create_random_planet(&planet);
		star = pdata->star;
		sdata = pdata->sdata;
	}
	void create_random_system(StellarSystem** system) {
/*		double min_age = sdata->evolution_ages[0];
		double precision = 0.01;
		double max_age = star->get_lifetime();
		OrbitSolver orb(min_age, max_age, precision,
				stellar_system_diff_eq, stellar_system_jacobian);*/
		*system = new StellarSystem(star, planet);
	}
};

///\brief Returns an evolutionary track for a quantity that is polynomial in 
///both mass and age.
PolynomialEvolutionTrack *exact_track(
		///The polynomial coefficients giving the mass and age dependence.
		///The first index should be age and the second one mass.
		const std::valarray< std::valarray<double> > &poly_coef, 
		
		///The stellar mass for which the tracks is necessary.
		double mass,
		
		///The power to scale ages for low mass stars.
		double low_mass_age_scaling=0,

		///The power to scale ages for high mass stars.
		double high_mass_age_scaling=0,
		
		///The mass on which to base the scaling.
		double scale_mass=NaN);

///\brief A stellar evolution class that initializes itself from a set of
///coefficients.
///
///\ingroup UnitTests_group
class PolynomialStellarEvolution : public StellarEvolution {
public:
	PolynomialStellarEvolution(
			///A list of the masses for which tracks will be evaluated and
			///passed to the StellarEvolution interpolation.
			const std::valarray<double> &masses,

			///A list of the masses for which tracks will be evaluated and
			///passed to the StellarEvolution interpolation.
			const std::valarray<double> &ages,

			///The coefficients of the stellar radius polynomial. First index
			///is age, the second is mass.
			const std::valarray< std::valarray<double> > &r_coef,

			///The coefficients of the convective zone moment of inertia 
			///polynomial. First index is age, the second is mass.
			const std::valarray< std::valarray<double> > &Iconv_coef,

			///The coefficients of the entire star moment of inertia 
			///polynomial. First index is age, the second is mass.
			const std::valarray< std::valarray<double> > &Itot_coef,

			///The coefficients of the mass of the radiative core polynomial.
			///First index is age, the second is mass.
			const std::valarray< std::valarray<double> > &Mrad_coef,

			///The coefficients of the radius of the radiative core
			///polynomial. First index is age, the second is mass.
			const std::valarray< std::valarray<double> > &Rcore_coef,

			///The power to scale ages for low mass stars.
			double low_mass_age_scaling,
			
			///The power to scale ages for high mass stars.
			double high_mass_age_scaling);
};

///\brief Represents a function of the form offset + scale*exp(rate*x) as
///well as its derivative.
///
///\ingroup UnitTests_group
class ExponentialPlusFunc : public OneArgumentDiffFunction,
	FunctionDerivatives {
private:
	double __scale, __rate, __deriv_x;
	OneArgumentDiffFunction *__offset;
public:
	///Create the function with the given parameters.
	ExponentialPlusFunc(OneArgumentDiffFunction *offset, double scale, double rate,
			double deriv_x=NaN) :
		__scale(scale), __rate(rate), __deriv_x(deriv_x), __offset(offset) {}

	///Evaluate the function at the given x.
	double operator()(double x) const
	{return (*__offset)(x) + __scale*std::exp(__rate*x);}

	///Returns the derivatives at the given x.
	const FunctionDerivatives *deriv(double x) const
	{return new ExponentialPlusFunc(__offset, __scale, __rate, x);}

	///\brief Return the given order-th derivative at x set by previous call
	///to deriv().
	double order(unsigned deriv_order=1) const;

	///The upper end of the range over which the function is defined.
	double range_high() const {return __offset->range_high();}

	///The lower end of the range over which the function is defined.
	double range_low() const {return __offset->range_low();}

	///Iterator over the points at which the function takes the given value.
	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};
///\brief Represents the sum of two functions and the derivative.
///
///\ingroup UnitTests_group
class FuncPlusFunc : public OneArgumentDiffFunction, FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f1, *__f2;
	double __deriv_x;
public:
	///Creates the function.
	FuncPlusFunc(const OneArgumentDiffFunction *f1,
			const OneArgumentDiffFunction *f2, double deriv_x=NaN) :
		__f1(f1), __f2(f2), __deriv_x(deriv_x) {}

	///Evaluates the function at the given x.
	double operator()(double x) const {return (*__f1)(x)+(*__f2)(x);}

	///Returns the derivative.
	const FunctionDerivatives *deriv(double x) const
	{return new FuncPlusFunc(__f1, __f2, x);}

	///For a derivative returns the given order.
	double order(unsigned deriv_order=1) const;

	///The upper end of the range where the function is defined.
	double range_high() const
	{return std::min(__f1->range_high(), __f2->range_high());}

	///The lower end of the range where the function is defined.
	double range_low() const
	{return std::max(__f1->range_low(), __f2->range_low());}

	///Iterator over the x values where the function takes the given value.
	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///\brief Several functions stiched together.
///
///\ingroup UnitTests_group
class PiecewiseFunction : public OneArgumentDiffFunction,FunctionDerivatives{
private:
	double __deriv_x, __range_low, __range_high;
	std::list<const OneArgumentDiffFunction *> __pieces;
public:
	///Create the function.
	PiecewiseFunction(const std::list<const OneArgumentDiffFunction *>
			&pieces=std::list<const OneArgumentDiffFunction *>(),
			double deriv_x=NaN);

	///\brief Adds another piece of the function, where it applies is defined
	///by its range.
	void add_piece(const OneArgumentDiffFunction *piece);

	///Evaluates the function at the given x.
	double operator()(double x) const;

	///Returns the derivatives at the given x.
	const FunctionDerivatives *deriv(double x) const
	{return new PiecewiseFunction(__pieces, x);}

	///For a derivative object returns the derivative of the given order.
	double order(unsigned deriv_order=1) const;

	///The upper end of the range over which the function is defined.
	double range_high() const {return __range_high;}

	///The lower end of the range over which the function is defined.
	double range_low() const {return __range_low;}

	///\brief An iterator over the x values where the function takes the
	///given value.
	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///\brief The ratio of two functions;
///
///\ingroup UnitTests_group
class FunctionRatio: public OneArgumentDiffFunction,FunctionDerivatives {
private:
	double __deriv_x;
	const OneArgumentDiffFunction *__f1, *__f2;
public:
	///Create the function
	FunctionRatio(const OneArgumentDiffFunction *f1,
			const OneArgumentDiffFunction *f2, double deriv_x=NaN) :
		__deriv_x(deriv_x), __f1(f1), __f2(f2) {}

	///Evaluates the function at the given x.
	double operator()(double x) const {return (*__f1)(x)/(*__f2)(x);}

	///Returns the derivatives at the given x.
	const FunctionDerivatives *deriv(double x) const
	{return new FunctionRatio(__f1, __f2, x);}

	///For a derivative object returns the derivative of the given order.
	double order(unsigned deriv_order=1) const;

	///The upper end of the range over which the function is defined.
	double range_high() const
	{return std::min(__f1->range_high(), __f2->range_high());}

	///The lower end of the range over which the function is defined.
	double range_low() const
	{return std::max(__f1->range_low(), __f2->range_low());}

	///\brief An iterator over the x values where the function takes the
	///given value.
	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///\brief A function raised to some power
///
///\ingroup UnitTests_group
class FunctionToPower : public OneArgumentDiffFunction, FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __power, __deriv_x;
public:
	///Create the function.
	FunctionToPower(const OneArgumentDiffFunction *f, double power,
			double deriv_x=NaN) :
		__f(f), __power(power), __deriv_x(deriv_x) {}

	///Evaluates the function at the given x.
	double operator()(double x) const
	{return std::pow((*__f)(x), __power);}

	///Returns the derivatives at the given x.
	const FunctionDerivatives *deriv(double x) const
	{return new FunctionToPower(__f, __power, x);}

	///For a derivative object returns the derivative of the given order.
	double order(unsigned deriv_order=1) const;

	///The upper end of the range over which the function is defined.
	double range_high() const {return __f->range_high();}

	///The lower end of the range over which the function is defined.
	double range_low() const {return __f->range_low();}

	///\brief An iterator over the x values where the function takes the
	///given value.
	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///\brief A function scaled by some constant
///
///\ingroup UnitTests_group
class ScaledFunction : public OneArgumentDiffFunction, FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __scale, __deriv_x;
public:
	///Create the function.
	ScaledFunction(const OneArgumentDiffFunction *f, double scale,
			double deriv_x=NaN) :
		__f(f), __scale(scale), __deriv_x(deriv_x) {}

	///Evaluates the function at the given x.
	double operator()(double x) const
	{return __scale*(*__f)(x);}

	///Returns the derivatives at the given x.
	const FunctionDerivatives *deriv(double x) const
	{return new ScaledFunction(__f, __scale, x);}

	///For a derivative object returns the derivative of the given order.
	double order(unsigned deriv_order=1) const;

	///The upper end of the range over which the function is defined.
	double range_high() const {return __f->range_high();}

	///The lower end of the range over which the function is defined.
	double range_low() const {return __f->range_low();}

	///\brief An iterator over the x values where the function takes the
	///given value.
	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///\brief The natural logarithm of a function.
///
///\ingroup UnitTests_group
class LogFunction : public OneArgumentDiffFunction, FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __deriv_x;
public: 
	///Create the function.
	LogFunction(const OneArgumentDiffFunction *f, double deriv_x=NaN) :
		__f(f), __deriv_x(deriv_x) {}

	///Evaluates the function at the given x.
	double operator()(double x) const {return std::log((*__f)(x));}

	///Returns the derivatives at the given x.
	const FunctionDerivatives *deriv(double x) const
	{return new LogFunction(__f, x);}

	///For a derivative object returns the derivative of the given order.
	double order(unsigned deriv_order=1) const;

	///The upper end of the range over which the function is defined.
	double range_high() const {return __f->range_high();}

	///The lower end of the range over which the function is defined.
	double range_low() const {return __f->range_low();}

	///\brief An iterator over the x values where the function takes the
	///given value.
	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///Fills the given valarray with a random set of polynomial coefficients.
void rand_poly_coef(std::valarray< std::valarray<double> > &poly_coef,
		double max_mass=-1);

///Returns a random set of polynomial coefficients
std::valarray< std::valarray<double> > rand_poly_coef(double max_mass=-1);

///Returns an array of the values of the track at the given ages.
std::valarray<double> tabulate_track(PolynomialEvolutionTrack *track,
		std::valarray<double> ages, unsigned deriv_order=0);

///\brief Returns the value of the polynomial with the given coefficients at
///the given mass and age.
double eval_poly(const std::valarray< std::valarray<double> > &poly_coef,
		double mass, double age, double low_mass_age_scaling=0,
		double high_mass_age_scaling=0, double scale_mass=NaN);

///\brief Outputs the mass and age polynomial defined by the given polynomial
///coefficients array
std::ostream &operator<<(std::ostream &os, 
		const std::valarray< std::valarray<double> > &poly_coef);

///A uniform random real value in the given range.
double rand_value(double min, double max);

///A uniform integer value in the given range.
int rand_value(int min, int max);

///Given n, m and (n)C(m) returns (n)C(m+1)
unsigned next_binom_coef(unsigned n, unsigned m, unsigned nCm);

///\brief Returns new polynomial coefficienst such that output
///polynomial(mass, age+age_offset)=input polynomial(mass, age)
std::valarray< std::valarray<double> > offset_age(
		const std::valarray< std::valarray<double> > &poly_coef,
		double age_offset);

///\brief Returns new polynomial coefficienst such that output
///polynomial(age+age_offset)=input polynomial(age)
std::valarray< std::valarray<double> > offset_age(
		const std::valarray< std::valarray<double> > &poly_coef,
		double age_offset);

///Solves f(x)=0 for x.
double solve(double guess_x, double abs_precision, double rel_precision,
		double (*f)(double x, void *params),
		double (*df) (double x, void *params),
		void (*fdf) (double x, void *params, double *f, double *df),
		void *params);

#endif
