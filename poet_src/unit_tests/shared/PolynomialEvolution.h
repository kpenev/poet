#ifndef __POLYNOMIAL_EVOLUTION_H
#define __POLYNOMIAL_EVOLUTION_H

#include "../Core/Functions.h"
#include "../StellarEvolution/EvolvingStellarQuantity.h"
#include "../StellarEvolution/Interpolator.h"

#include "Common.h"
//#include "../StellarSystem.h"
//#include "../Star.h"
//#include "../Planet.h"
//#include "../OrbitSolver.h"
#include <valarray>
#include <gsl/gsl_roots.h>
//#include "YRECIO.h"
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace StellarEvolution {

    ///\brief An EvolvingStellar quantity that uses a polynomial instead of
    ///interpolating.
    ///
    ///It also serves as its own derivative.
    ///
    ///\ingroup UnitTests_group
    class PolynomialEvolutionQuantity : public EvolvingStellarQuantity,
        Core::FunctionDerivatives {
    private:
        ///The coefficients of the polynomial giving the evolution.
        std::valarray<double> __poly_coef;

        double 
            ///\brief The lower limit of the age at which this quantity can
            ///be evaluated.
            __xmin,

            ///\brief The upper limit of the age at which this quantity
            ///can be evaluated.
            __xmax;

        ///The location at which the derivative has been requested.
        mutable	double __deriv_x;

        ///An empty vector to serve as the list of discontinuities.
        std::vector<double> __empty_vector;

    public:
        ///Create an evolving quantity or its derivative.
        PolynomialEvolutionQuantity(
            ///The coefficients defining the polynomial that gives the
            ///evolution of the quantity.
            const std::valarray<double> &coefficients,

            ///The lower end of the range for which the quantity is defined.
            double range_low,

            ///The upper end of the range for which the quantity is defined.
            double range_high,

            ///The abscissa at which the derivatives are to be evaluated.
            double derivative_x=NaN
        ) :
            __poly_coef(coefficients),
            __xmin(range_low),
            __xmax(range_high),
            __deriv_x(derivative_x)
            {}

        void select_interpolation_region(double) const {}

        ///\brief Evaluates the polynomial at the given age.
        ///
        ///Throws an exception if age is outside the range over which the
        ///quantity is defined.
        double operator()(double x) const {__deriv_x = x; return order(0);}

        ///The lower end of the range where the function is defined
        double range_low() const {return __xmin;}

        ///The upper end of the range where the function is defined
        double range_high() const {return __xmax;}

        ///The ages at which the quantity may be discontinuous.
        virtual const std::vector<double> &discontinuities() const
        {return __empty_vector;}

        ///\brief See
        ///StellarEvolution::EvolvingStellarQuantity::prevous_discontinuity()
        virtual double previous_discontinuity() const {return __xmin;}

        ///\brief See
        ///StellarEvolution::EvolvingStellarQuantity::next_discontinuity()
        virtual double next_discontinuity() const {return __xmax;}

        ///\brief Set up the interpolation over the next interpolation region
        ///(between consecutive discontinuities.)
        virtual void enable_next_interpolation_region() const
        {assert(false);}

        ///The derivatives.
        const Core::FunctionDerivatives *deriv(double x) const
        {return new PolynomialEvolutionQuantity(__poly_coef,
                                                __xmin,
                                                __xmax,
                                                x);}

        ///The order-th derivative
        double order(unsigned deriv_order = 1) const;

        ///Prints an expression of the polynomial that this track represents.
        friend std::ostream &operator<<(
            std::ostream &os,
            const PolynomialEvolutionQuantity &track
        );

        ///A string identifying the type of quantity this is.
        std::string kind() const {return "PolynomialEvolutionQuantity";}
    };//End PolynomialEvolutionQuantity class.

    ///\brief Implements a StellarEvolution based on polynomial evolution
    ///quantities.
    ///
    ///\ingroup UnitTests_group
    class MockStellarEvolution : public StellarEvolution::Interpolator {
    private:
        std::valarray< std::valarray<double> >
            ///\brief The polynomial coefficients of the radius dependence on
            ///mass and age.
            __R,

            ///\brief The polynomial coefficients of the convective zone
            ///moment of inertia dependence on stellar mass and age.
            __Iconv,

            ///\brief The polynomial coefficients of the radiative zone
            ///moment of inertia dependence on stellar mass and age.
            __Irad,

            ///\brief The polynomial coefficients of total stellar moment of
            ///inertia dependence on stellar mass and age.
            __Itot,

            ///\brief The polynomial coefficients of the core radius
            ///dependence on stellar mass and age.
            __Rcore,

            ///\brief The polynomial coefficients of the core mass dependence
            ///on stellar mass and age.
            __Mcore,

            ///\brief The polynomial coefficients of the envelope mass
            ///dependence on stellar mass and age.
            __Menv,
            
            ///\brief The polynomial coefficients of the luminosity
            ///dependence on stellar mass and age.
            __Lum;

        ///The age in Gyr at which the core forms.
        double __core_formation_age;
    public:
        ///\brief Uses the given values and polynomials for the quantities.
        ///
        ///If empty polynomial coefficients are specified or NaN for constant
        ///quantities in which case random polynomials or values are used.
        MockStellarEvolution(
            ///The age at which the core forms in Gyr.
            double core_formation_age = Core::NaN,

            ///The polynomial coefficients of the radius dependence on mass
            ///and age.
            const std::valarray< std::valarray<double> > &
                R = std::valarray< std::valarray<double> >(),

            ///The polynomial coefficients of the convective zone moment of
            ///inertia dependence on stellar mass and age.
            const std::valarray< std::valarray<double> > &
                Iconv = std::valarray< std::valarray<double> >(),

            ///The polynomial coefficients of the radiative zone moment of
            ///inertia dependence on stellar mass and age.
            const std::valarray< std::valarray<double> > &
                Irad = std::valarray< std::valarray<double> >(),

            ///The polynomial coefficients of the core radius dependence on
            ///stellar mass and age.
            const std::valarray< std::valarray<double> > &
                Rcore = std::valarray< std::valarray<double> >(),

            ///The polynomial coefficients of the core mass dependence on
            ///stellar mass and age.
            const std::valarray< std::valarray<double> > &
                Mcore = std::valarray< std::valarray<double> >(),

            ///The polynomial coefficients of the luminosity dependence on
            ///stellar mass and age.
            const std::valarray< std::valarray<double> > &
                Lum = std::valarray< std::valarray<double> >()
        );

        ///See StellarEvolutino::Interpolator::operator()()
        EvolvingStellarQuantity *operator()(QuantityID quantity,
                                            double mass,
                                            double feh) const;

        ///The age in Gyr at which the core forms.
        double core_formation_age() const {return __core_formation_age;}
    };

    MockStellarEvolution *make_no_evolution(double Rstar = 1.0,
                                            double Iconv = 1.0);

    MockStellarEvolution *make_linear_I_evolution();

}//End StellarEvolution namespace.

#if 0
    ///\brief A messy class able to create random stars.
    ///
    ///\ingroup UnitTests_group
    class StarData {
    private:
        unsigned __num_stars_created;
    public:
        double mass,
               feh,
               radius,
               age,
               conv_spin,
               rad_spin,
               tidal_Q,
               wind_strength,
               wind_sat_freq,
               coupling_timescale,
               Q_trans_width,
               disk_lock_w,
               disk_lock_time;

        PolynomialEvolutionTrack *Lrad_track,
                                 *Lconv_track;

        std::valarray<double> evolution_masses,
                              evolution_ages,
                              Lrad,
                              Lconv,
                              Iconv,
                              Irad,
                              Mrad_deriv,
                              Lconv_deriv,
                              Lrad_deriv,
                              Rrad,
                              all_radii;

        std::valarray< std::valarray<double> > r_coef,
                                               Iconv_coef,
                                               Irad_coef,
                                               Mrad_coef,
                                               Rcore_coef;

        StarData();

        ///Creates a star with random properties (mass, radius, age, ...,
        ///random polynomial time functions as moments of inertia ...).
        void create_random_star(Star::EvolvingStar **star);

        ///Deletes the convective and radiative angular momenta tracks if
        ///they were allocated.
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
		*system = new StellarSystem(*star, *planet, "");
	}
};
#endif

namespace StellarEvolution {

    ///\brief Returns an evolutionary track for a quantity that is polynomial
    ///in both mass and age.
    PolynomialEvolutionQuantity *exact_track(
        ///The polynomial coefficients giving the mass and age dependence.
        ///The first index should be age and the second one mass.
        const std::valarray< std::valarray<double> > &poly_coef, 

        ///The stellar mass for which the tracks is necessary.
        double mass,

        ///The power to scale ages for low mass stars.
        double low_mass_age_scaling = 0,

        ///The power to scale ages for high mass stars.
        double high_mass_age_scaling = 0,

        ///The mass on which to base the scaling.
        double scale_mass = NaN
    );

} //End StellarEvolution namespace.

#if 0
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
#endif

///\brief Represents a function of the form offset + scale*exp(rate*x) as
///well as its derivative.
///
///\ingroup UnitTests_group
class ExponentialPlusFunc : public Core::OneArgumentDiffFunction,
	                               Core::FunctionDerivatives
{
private:
	double __scale, __rate, __deriv_x;
    Core::OneArgumentDiffFunction *__offset;
public:
	///Create the function with the given parameters.
	ExponentialPlusFunc(Core::OneArgumentDiffFunction *offset,
                        double scale,
                        double rate,
                        double deriv_x = Core::NaN) :
		__scale(scale),
        __rate(rate),
        __deriv_x(deriv_x),
        __offset(offset)
    {}

	///Evaluate the function at the given x.
	double operator()(double x) const
	{return (*__offset)(x) + __scale * std::exp(__rate * x);}

	///Returns the derivatives at the given x.
	const Core::FunctionDerivatives *deriv(double x) const
	{return new ExponentialPlusFunc(__offset, __scale, __rate, x);}

	///\brief Return the given order-th derivative at x set by previous call
	///to deriv().
	double order(unsigned deriv_order = 1) const;

	///The upper end of the range over which the function is defined.
	double range_high() const {return __offset->range_high();}

	///The lower end of the range over which the function is defined.
	double range_low() const {return __offset->range_low();}

	///Iterator over the points at which the function takes the given value.
    Core::InterpSolutionIterator crossings(double) const
	{return Core::InterpSolutionIterator();}
};

///\brief Represents the sum of two functions and the derivative.
///
///\ingroup UnitTests_group
class FuncPlusFunc : public Core::OneArgumentDiffFunction,
                            Core::FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f1, *__f2;
	double __deriv_x;
public:
	///Creates the function.
	FuncPlusFunc(const OneArgumentDiffFunction *f1,
                 const OneArgumentDiffFunction *f2,
                 double deriv_x = Core::NaN) :
		__f1(f1), __f2(f2), __deriv_x(deriv_x)
    {}

	///Evaluates the function at the given x.
	double operator()(double x) const {return (*__f1)(x) + (*__f2)(x);}

	///Returns the derivative.
	const FunctionDerivatives *deriv(double x) const
	{return new FuncPlusFunc(__f1, __f2, x);}

	///For a derivative returns the given order.
	double order(unsigned deriv_order = 1) const;

	///The upper end of the range where the function is defined.
	double range_high() const
	{return std::min(__f1->range_high(), __f2->range_high());}

	///The lower end of the range where the function is defined.
	double range_low() const
	{return std::max(__f1->range_low(), __f2->range_low());}

	///Iterator over the x values where the function takes the given value.
    Core::InterpSolutionIterator crossings(double) const
	{return Core::InterpSolutionIterator();}
};

///\brief Several functions stiched together.
///
///\ingroup UnitTests_group
class PiecewiseFunction : public Core::OneArgumentDiffFunction,
                                 Core::FunctionDerivatives{
private:
	double __deriv_x, __range_low, __range_high;
	std::list<const OneArgumentDiffFunction *> __pieces;
public:
	///Create the function.
	PiecewiseFunction(const std::list<const OneArgumentDiffFunction *>
                      &pieces = std::list<const OneArgumentDiffFunction *>(),
                      double deriv_x = Core::NaN);

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
    Core::InterpSolutionIterator crossings(double) const
	{return Core::InterpSolutionIterator();}
};

///\brief The ratio of two functions;
///
///\ingroup UnitTests_group
class FunctionRatio: public Core::OneArgumentDiffFunction,
                            Core::FunctionDerivatives {
private:
	double __deriv_x;
	const OneArgumentDiffFunction *__f1, *__f2;
public:
	///Create the function
	FunctionRatio(const OneArgumentDiffFunction *f1,
                  const OneArgumentDiffFunction *f2,
                  double deriv_x = Core::NaN) :
		__deriv_x(deriv_x), __f1(f1), __f2(f2) {}

	///Evaluates the function at the given x.
	double operator()(double x) const {return (*__f1)(x) / (*__f2)(x);}

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
    Core::InterpSolutionIterator crossings(double) const
	{return Core::InterpSolutionIterator();}
};

///\brief A function raised to some power
///
///\ingroup UnitTests_group
class FunctionToPower : public Core::OneArgumentDiffFunction,
                               Core::FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __power, __deriv_x;
public:
	///Create the function.
	FunctionToPower(const OneArgumentDiffFunction *f,
                    double power,
                    double deriv_x = Core::NaN) :
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
    Core::InterpSolutionIterator crossings(double) const
	{return Core::InterpSolutionIterator();}
};

///\brief A function scaled by some constant
///
///\ingroup UnitTests_group
class ScaledFunction : public Core::OneArgumentDiffFunction,
                              Core::FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __scale, __deriv_x;
public:
	///Create the function.
	ScaledFunction(const OneArgumentDiffFunction *f,
                   double scale,
                   double deriv_x = Core::NaN) :
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
    Core::InterpSolutionIterator crossings(double) const
	{return Core::InterpSolutionIterator();}
};

///\brief The natural logarithm of a function.
///
///\ingroup UnitTests_group
class LogFunction : public Core::OneArgumentDiffFunction,
                           Core::FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __deriv_x;
public: 
	///Create the function.
	LogFunction(const OneArgumentDiffFunction *f,
                double deriv_x = Core::NaN) :
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
    Core::InterpSolutionIterator crossings(double) const
	{return Core::InterpSolutionIterator();}
};

///\brief The cosine of a function.
///
///\ingroup UnitTests_group
class CosFunction : public Core::OneArgumentDiffFunction,
                           Core::FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __deriv_x;
public: 
	///Create the function.
	CosFunction(const OneArgumentDiffFunction *f,
                double deriv_x = Core::NaN) :
		__f(f), __deriv_x(deriv_x) {}

	///Evaluates the function at the given x.
	double operator()(double x) const {return std::cos((*__f)(x));}

	///Returns the derivatives at the given x.
	const FunctionDerivatives *deriv(double x) const
	{return new CosFunction(__f, x);}

	///For a derivative object returns the derivative of the given order.
	double order(unsigned deriv_order=1) const;

	///The upper end of the range over which the function is defined.
	double range_high() const {return __f->range_high();}

	///The lower end of the range over which the function is defined.
	double range_low() const {return __f->range_low();}

	///\brief An iterator over the x values where the function takes the
	///given value.
    Core::InterpSolutionIterator crossings(double) const
	{return Core::InterpSolutionIterator();}
};


#if 0
///Returns an array of the values of the track at the given ages.
std::valarray<double> tabulate_track(PolynomialEvolutionTrack *track,
		std::valarray<double> ages, unsigned deriv_order=0);

///\brief Returns the value of the polynomial with the given coefficients at
///the given mass and age.
double eval_poly(const std::valarray< std::valarray<double> > &poly_coef,
                 double mass, double age, double low_mass_age_scaling=0,
		double high_mass_age_scaling=0, double scale_mass=NaN);
#endif

///Solves f(x)=0 for x.
double solve(double guess_x,
             double abs_precision,
             double rel_precision,
             double (*f)(double x, void *params),
             double (*df) (double x, void *params),
             void (*fdf) (double x, void *params, double *f, double *df),
             void *params);

#endif
