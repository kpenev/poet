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

///Returns the fractional difference between x and y.
bool check_diff(double x, double y, double frac_tolerance, 
		double abs_tolerance);

bool check_diff(std::valarray<double> x, std::valarray<double> y,
		std::valarray<double> frac_tolerance,
		std::valarray<double> abs_tolerance);

bool isEqual(double a, double b);
double getError(double predicted, double actual);
bool approxEqual(double predicted, double actual, double thres=0.02);

const double min_stellar_mass=0.4, max_low_mass=1.075, max_stellar_mass=1.3,
	  min_age=1e-7, max_age=10.0;
const double min_planet_mass=10, max_planet_mass=80, min_planet_radius=5,
		max_planet_radius=15;

std::ostream &operator<<(std::ostream &os,
		const std::valarray<double> &array);


class PolynomialEvolutionTrack : public EvolvingStellarQuantity,
	FunctionDerivatives {
private:
	std::valarray<double> poly_coef;
	double xmin, xmax;
	mutable	double deriv_x;
public:
	PolynomialEvolutionTrack(const std::valarray<double> &coefficients,
			double range_low, double range_high, double derivative_x=NaN) :
		poly_coef(coefficients), xmin(range_low), xmax(range_high),
		deriv_x(derivative_x) {}

	///Evaluates the polynomial at the given x if it is in the range,
	///throws an exception otherwise.
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

	std::string kind() const {return "PolynomialEvolutionTrack";}
};

class MockStellarEvolution : public StellarEvolution {
private:
	std::valarray< std::valarray<double> > R, Iconv, Irad, Itot, Rcore,
		Mcore, Menv;

	double core_formation;
public:
	///Uses the given values and polynomials for the quantities, unless they
	///are empty or NaN in which case random polynomials or values are used.
	MockStellarEvolution(double core_formation_age=NaN,
			const std::valarray< std::valarray<double> > &radius=
				std::valarray< std::valarray<double> >(),
			const std::valarray< std::valarray<double> >
				&conv_moment_of_inertia=
				std::valarray< std::valarray<double> >(),
			const std::valarray< std::valarray<double> > 
				&rad_moment_of_inertia=
				std::valarray< std::valarray<double> >(),
			const std::valarray< std::valarray<double> > &core_radius=
				std::valarray< std::valarray<double> >(),
			const std::valarray< std::valarray<double> > &core_mass=
				std::valarray< std::valarray<double> >());

	///Returns a single argument function which gives the moment of
	///inertia of the specified zone of a star of the specified mass as a
	///function of age. The result must be destroyed when it becomes
	///obsolete. If the present age and radius are specified, the result is
	///scaled by (present_radius/r(present_age))^2.
	const EvolvingStellarQuantity *interpolate_moment_of_inertia(
			double stellar_mass, StellarZone zone, double present_age=-1)
		const;

	///Returns a single argument function which gives the radius of a 
	///star of the specified mass as a function of age. The result must
	///be destroyed when it becomes obsolete. If the present age and radius
	///are specified, the result is scaled by
	///(present_radius/r(present_age)).
	const EvolvingStellarQuantity *interpolate_radius(
			double stellar_mass, double present_age=-1) const;

	///Returns a single argument function which gives the mass of
	///of the specified zone of a star of the specified mass as a
	///function of age. The result must be destroyed when it becomes
	///obsolete.
	const EvolvingStellarQuantity *interpolate_zone_mass(
			double stellar_mass, StellarZone zone) const;

	///Returns a single argument function which gives the radius of the 
	///convective-radiative boundary for a star of the specified mass as
	///a function of age. The result must be destroyed when it becomes 
	///obsolete. If the present age and radius are specified, the result is
	///scaled by (present_radius/r(present_age)).
	const EvolvingStellarQuantity *interpolate_core_boundary(
			double stellar_mass, double present_age=-1) const;

	double core_formation_age() const {return core_formation;}
};


///A class able to create random stars.
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

///Returns an evolutionary track for a quantity that is polynomial in 
///both mass and age, with the given polynomial coefficients.
///The first index should be age and the second one mass.
PolynomialEvolutionTrack *exact_track(
		const std::valarray< std::valarray<double> > &poly_coef, 
		double mass, double low_mass_age_scaling=0,
		double high_mass_age_scaling=0, double scale_mass=NaN);

///A stellar evolution class that initializes itself from a set of
///coefficients
class PolynomialStellarEvolution : public StellarEvolution {
public:
	PolynomialStellarEvolution(
			const std::valarray<double> &masses,
			const std::valarray<double> &ages,
			const std::valarray< std::valarray<double> > &r_coef,
			const std::valarray< std::valarray<double> > &Iconv_coef,
			const std::valarray< std::valarray<double> > &Itot_coef,
			const std::valarray< std::valarray<double> > &Mrad_coef,
			const std::valarray< std::valarray<double> > &Rcore_coef,
			double low_mass_age_scaling, double high_mass_age_scaling);

};

///Represents the function offset + scale*exp(rate*x)
class ExponentialPlusFunc : public OneArgumentDiffFunction,
	FunctionDerivatives {
private:
	double __scale, __rate, __deriv_x;
	OneArgumentDiffFunction *__offset;
public:
	ExponentialPlusFunc(OneArgumentDiffFunction *offset, double scale, double rate,
			double deriv_x=NaN) :
		__scale(scale), __rate(rate), __deriv_x(deriv_x), __offset(offset) {}

	double operator()(double x) const
	{return (*__offset)(x) + __scale*std::exp(__rate*x);}

	const FunctionDerivatives *deriv(double x) const
	{return new ExponentialPlusFunc(__offset, __scale, __rate, x);}

	double order(unsigned deriv_order=1) const;

	double range_high() const {return __offset->range_high();}
	double range_low() const {return __offset->range_low();}

	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///Represents the sum of two functions
class FuncPlusFunc : public OneArgumentDiffFunction, FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f1, *__f2;
	double __deriv_x;
public:
	FuncPlusFunc(const OneArgumentDiffFunction *f1,
			const OneArgumentDiffFunction *f2, double deriv_x=NaN) :
		__f1(f1), __f2(f2), __deriv_x(deriv_x) {}

	double operator()(double x) const {return (*__f1)(x)+(*__f2)(x);}

	const FunctionDerivatives *deriv(double x) const
	{return new FuncPlusFunc(__f1, __f2, x);}

	double order(unsigned deriv_order=1) const;

	double range_high() const
	{return std::min(__f1->range_high(), __f2->range_high());}

	double range_low() const
	{return std::max(__f1->range_low(), __f2->range_low());}

	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

class PiecewiseFunction : public OneArgumentDiffFunction,FunctionDerivatives{
private:
	double __deriv_x, __range_low, __range_high;
	std::list<const OneArgumentDiffFunction *> __pieces;
public:
	PiecewiseFunction(const std::list<const OneArgumentDiffFunction *>
			&pieces=std::list<const OneArgumentDiffFunction *>(),
			double deriv_x=NaN);

	void add_piece(const OneArgumentDiffFunction *piece);

	double operator()(double x) const;

	const FunctionDerivatives *deriv(double x) const
	{return new PiecewiseFunction(__pieces, x);}

	double order(unsigned deriv_order=1) const;

	double range_high() const {return __range_high;}

	double range_low() const {return __range_low;}

	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///The ratio of two functions;
class FunctionRatio: public OneArgumentDiffFunction,FunctionDerivatives {
private:
	double __deriv_x;
	const OneArgumentDiffFunction *__f1, *__f2;
public:
	FunctionRatio(const OneArgumentDiffFunction *f1,
			const OneArgumentDiffFunction *f2, double deriv_x=NaN) :
		__deriv_x(deriv_x), __f1(f1), __f2(f2) {}

	double operator()(double x) const {return (*__f1)(x)/(*__f2)(x);}

	const FunctionDerivatives *deriv(double x) const
	{return new FunctionRatio(__f1, __f2, x);}

	double order(unsigned deriv_order=1) const;

	double range_high() const
	{return std::min(__f1->range_high(), __f2->range_high());}

	double range_low() const
	{return std::max(__f1->range_low(), __f2->range_low());}

	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///A function raised to some power
class FunctionToPower : public OneArgumentDiffFunction, FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __power, __deriv_x;
public:
	FunctionToPower(const OneArgumentDiffFunction *f, double power,
			double deriv_x=NaN) :
		__f(f), __power(power), __deriv_x(deriv_x) {}

	double operator()(double x) const
	{return std::pow((*__f)(x), __power);}

	const FunctionDerivatives *deriv(double x) const
	{return new FunctionToPower(__f, __power, x);}

	double order(unsigned deriv_order=1) const;

	double range_high() const {return __f->range_high();}

	double range_low() const {return __f->range_low();}

	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

///A function scaled by some constant
class ScaledFunction : public OneArgumentDiffFunction, FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __scale, __deriv_x;
public:
	ScaledFunction(const OneArgumentDiffFunction *f, double scale,
			double deriv_x=NaN) :
		__f(f), __scale(scale), __deriv_x(deriv_x) {}

	double operator()(double x) const
	{return __scale*(*__f)(x);}

	const FunctionDerivatives *deriv(double x) const
	{return new ScaledFunction(__f, __scale, x);}

	double order(unsigned deriv_order=1) const;

	double range_high() const {return __f->range_high();}

	double range_low() const {return __f->range_low();}

	InterpSolutionIterator crossings(double) const
	{return InterpSolutionIterator();}
};

class LogFunction : public OneArgumentDiffFunction, FunctionDerivatives {
private:
	const OneArgumentDiffFunction *__f;
	double __deriv_x;
public: 
	LogFunction(const OneArgumentDiffFunction *f, double deriv_x=NaN) :
		__f(f), __deriv_x(deriv_x) {}

	double operator()(double x) const {return std::log((*__f)(x));}

	const FunctionDerivatives *deriv(double x) const
	{return new LogFunction(__f, x);}

	double order(unsigned deriv_order=1) const;

	double range_high() const {return __f->range_high();}

	double range_low() const {return __f->range_low();}

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

///Returns the value of the polynomial with the given coefficients at the
///given mass and age.
double eval_poly(const std::valarray< std::valarray<double> > &poly_coef,
		double mass, double age, double low_mass_age_scaling=0,
		double high_mass_age_scaling=0, double scale_mass=NaN);

///Outputs the mass and age polynomial defined by the given polynomial
///coefficients array
std::ostream &operator<<(std::ostream &os, 
		const std::valarray< std::valarray<double> > &poly_coef);

double rand_value(double min, double max);

int rand_value(int min, int max);

///Given n, m and (n)C(m) returns (n)C(m+1)
unsigned next_binom_coef(unsigned n, unsigned m, unsigned nCm);

///Returns new polynomial coefficienst such that
///output polynomial(mass, age+age_offset)=input polynomial(mass, age)
std::valarray< std::valarray<double> > offset_age(
		const std::valarray< std::valarray<double> > &poly_coef,
		double age_offset);

///Returns new polynomial coefficienst such that
///output polynomial(age+age_offset)=input polynomial(age)
std::valarray< std::valarray<double> > offset_age(
		const std::valarray< std::valarray<double> > &poly_coef,
		double age_offset);

///Solves the given equation and its derivative using
double solve(double guess_x, double abs_precision, double rel_precision,
		double (*f)(double x, void *params),
		double (*df) (double x, void *params),
		void (*fdf) (double x, void *params, double *f, double *df),
		void *params);

#endif
