#include "Common.h"
#include "YRECIO.h"
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
///Returns the fractional difference between x and y.
bool check_diff(double x, double y, double frac_tolerance, 
		double abs_tolerance)
{
	if(std::isnan(x)) return std::isnan(y);
	return std::abs(x-y)<=(abs_tolerance+
		frac_tolerance*std::max(std::abs(x), std::abs(y)));
}

bool check_diff(std::valarray<double> x, std::valarray<double> y,
		std::valarray<double> frac_tolerance,
		std::valarray<double> abs_tolerance)
{
	bool result=true;
	for(size_t i=0; i<x.size(); i++)
		result=result && (std::abs(x[i]-y[i])<=
				(abs_tolerance[i]+
				frac_tolerance[i]*std::max(std::abs(x[i]), std::abs(y[i]))));
	return result;
}

std::ostream &operator<<(std::ostream &os,
		const std::valarray<double> &array)
{
	os << array[0];
	for(size_t i=1; i<array.size(); i++)
		os << ", " << array[i];
	return os;
}

bool isEqual(double a, double b) {
	return std::abs(a-b) < 1e-20;
}

double getError(double predicted, double actual) {
	return std::abs((actual - predicted)/actual);
}

bool approxEqual(double predicted, double actual, double thres) {
	if (isEqual(predicted, actual)) return true;
	return getError(predicted, actual) < 0.02;
}

///The order-th derivative
double PolynomialEvolutionTrack::order(unsigned deriv_order) const
{
	if(deriv_x<xmin || deriv_x>xmax) {
		std::ostringstream msg;
		msg << "Attempt to evaluate a PolynomialEvolutionTrack defined over "
			<< xmin << " < x < " << xmax << " at x=" << deriv_x
			<< ", polynomial coefficients:" << poly_coef;
		std::cout << msg.str() << std::endl;
		*static_cast<int*>(NULL)=0;
		throw Error::Runtime(msg.str());
	}
	double xn=1, result=0;
	for(size_t n=deriv_order; n<poly_coef.size(); n++) {
		double coef_mod=1;
		for(size_t i=0; i<deriv_order; i++) coef_mod*=(n-i);
		result+=coef_mod*poly_coef[n]*xn;
		xn*=deriv_x;
	}
	return result;
}

///Uses the given values and polynomials for the quantities, unless they
///are empty or NaN in which case random polynomials or values are used.
MockStellarEvolution::MockStellarEvolution(double core_formation_age,
		const std::valarray< std::valarray<double> > &radius,
		const std::valarray< std::valarray<double> > &conv_moment_of_inertia,
		const std::valarray< std::valarray<double> > &rad_moment_of_inertia,
		const std::valarray< std::valarray<double> > &core_radius,
		const std::valarray< std::valarray<double> > &core_mass) :
	R(radius), Iconv(conv_moment_of_inertia), Irad(rad_moment_of_inertia),
	Rcore(core_radius), Mcore(core_mass), core_formation(core_formation_age)
{
	if(std::isnan(core_formation))
		core_formation=rand_value(min_age, 2.1);
	if(R.size()==0) rand_poly_coef(R);
	if(Iconv.size()==0) rand_poly_coef(Iconv);
	if(Irad.size()==0) rand_poly_coef(Irad);
	if(Rcore.size()==0) rand_poly_coef(Rcore);
	if(Mcore.size()==0) rand_poly_coef(Mcore);
	if(core_formation_age>0) {
		Irad[0]=0;
		Rcore[0]=0;
		Mcore[0]=0;
		Irad=offset_age(Irad, core_formation);
		Rcore=offset_age(Rcore, core_formation);
		Mcore=offset_age(Mcore, core_formation);
	}
	Menv.resize(Mcore.size(), std::valarray<double>(Mcore[0].size())); 
	Itot.resize(Iconv.size(), std::valarray<double>(Iconv[0].size())); 
	Itot=Iconv+Irad;
	Menv-=Mcore;
	Menv[0][1]+=1;
}

///Returns a single argument function which gives the moment of
///inertia of the specified zone of a star of the specified mass as a
///function of age. The result must be destroyed when it becomes
///obsolete. If the present age and radius are specified, the result is
///scaled by (present_radius/r(present_age))^2.
const EvolvingStellarQuantity *
MockStellarEvolution::interpolate_moment_of_inertia(
		double stellar_mass, StellarZone zone, double present_age) const
{
	switch(zone) {
		case convective : return exact_track(Iconv, stellar_mass);
		case radiative : return exact_track(Irad, stellar_mass);
		case total : return exact_track(Itot, stellar_mass);
		default : throw Error::BadStellarZone("Invalid stellar zone in "
				  "MockStellarEvolution::interpolate_moment_of_inertia.");
	}
}

///Returns a single argument function which gives the radius of a 
///star of the specified mass as a function of age. The result must
///be destroyed when it becomes obsolete. If the present age and radius
///are specified, the result is scaled by
///(present_radius/r(present_age)).
const EvolvingStellarQuantity *MockStellarEvolution::interpolate_radius(
		double stellar_mass, double present_age) const
{
	return exact_track(R, stellar_mass);
}

///Returns a single argument function which gives the mass of
///of the specified zone of a star of the specified mass as a
///function of age. The result must be destroyed when it becomes
///obsolete.
const EvolvingStellarQuantity *MockStellarEvolution::interpolate_zone_mass(
		double stellar_mass, StellarZone zone) const
{
	switch(zone) {
		case convective : return exact_track(Menv, stellar_mass);
		case radiative : return exact_track(Mcore, stellar_mass);
		default : throw Error::BadStellarZone("Only convectie and radiative "
						  "zone masses can be requested in "
						  "MockStellarEvolution::interpolate_zone_mass.");
	}
}

///Returns a single argument function which gives the radius of the 
///convective-radiative boundary for a star of the specified mass as
///a function of age. The result must be destroyed when it becomes 
///obsolete. If the present age and radius are specified, the result is
///scaled by (present_radius/r(present_age)).
const EvolvingStellarQuantity *
MockStellarEvolution::interpolate_core_boundary(
		double stellar_mass, double present_age) const
{
	return exact_track(Rcore, stellar_mass);
}

StarData::StarData() :
	num_stars_created(0), Lrad_track(NULL), Lconv_track(NULL),
	evolution_masses(100), evolution_ages(100), r_coef(rand_poly_coef()),
	Iconv_coef(rand_poly_coef()), Itot_coef(rand_poly_coef()),
	Mrad_coef(rand_poly_coef()), Rcore_coef(rand_poly_coef())
{
	for (size_t age_ind=0; age_ind < evolution_ages.size(); age_ind++) {
		evolution_ages[age_ind]=min_age+
			age_ind*(max_age-min_age)/evolution_ages.size();
	}

	for (size_t mass_ind=0; mass_ind < evolution_masses.size(); mass_ind++) {
		evolution_masses[mass_ind] = min_stellar_mass +
			mass_ind*(max_stellar_mass - min_stellar_mass)/
			evolution_masses.size();
	}
}

///Deletes the convective and radiative angular momenta tracks if they
///were allocated.
StarData::~StarData()
{
	if(Lrad_track) delete Lrad_track;
	if(Lconv_track) delete Lconv_track;
}

void StarData::create_random_star(Star **star)
{
	mass=rand_value(min_stellar_mass, max_stellar_mass);
	age=rand_value(min_age, 2.1);
	radius=(rand_value(0.9, 1.1))*eval_poly(r_coef, mass, age);
	conv_spin=rand_value(0.0, 1.0);
	rad_spin=rand_value(0.0, 1.0);
	tidal_Q=std::pow(10.0, rand_value(5.0,10.0));
	wind_strength=rand_value(0.0, 1.0);
	wind_sat_freq=rand_value(0.0, 1.0);
	disk_lock_w=rand_value(0.0, wind_sat_freq);
	disk_lock_time=evolution_ages[rand_value(0, evolution_ages.size()-1)];
	coupling_timescale = rand_value(0.0, 1.0);
	Q_trans_width = 0;

	//(*star) = new Star("test2", age, mass, radius, conv_spin, rad_spin, tidal_Q,
	//	wind_strength, wind_sat_freq, coupling_timescale,
	//Q_trans_width, PolynomialStellarEvolution(evolution_masses,
	//evolution_ages, r_coef, Iconv_coef, Itot_coef, Mrad_coef,
	//Rcore_coef,0,0));
#ifdef NO_SERIALIZE
//	YRECEvolution evol("../YREC",0, 2, 3);
//	evol.save_state("test_save");
	MockStellarEvolution evol;
	(*star) = new Star(mass, tidal_Q, wind_strength, wind_sat_freq,
			coupling_timescale, Q_trans_width, disk_lock_w, disk_lock_time,
			evol, age, conv_spin, rad_spin);
#endif
#ifndef NO_SERIALIZE
	YRECEvolution evol;
	evol.load_state("test_save");
	(*star) = new Star(mass, tidal_Q, wind_strength, wind_sat_freq,
			coupling_timescale, Q_trans_width, disk_lock_w, disk_lock_time,
			evol, age, conv_spin, rad_spin);
#endif

	//(*star) = new Star("test2", age, mass, radius, conv_spin, rad_spin, tidal_Q,
	//			wind_strength, wind_sat_freq, coupling_timescale,
	//		Q_trans_width, evol);
	//initialize Lrad, Lconv,Iconv, Irad, Mrad_deriv, Lconv_deriv,
	//Lrad_deriv, Rrad;
	Lconv_track=exact_track(rand_poly_coef(), mass);
	Lrad_track=exact_track(rand_poly_coef(), mass);
	Lconv.resize(evolution_ages.size());
	Lconv = tabulate_track(Lconv_track, evolution_ages);
	Lrad.resize(evolution_ages.size());
	Lrad = tabulate_track(Lrad_track, evolution_ages);
	Iconv.resize(evolution_ages.size());
	Iconv = tabulate_track(exact_track(Iconv_coef, mass), evolution_ages);
	std::valarray<double> Itot =
		tabulate_track(exact_track(Itot_coef, mass), evolution_ages);

	std::list<double> Irad_list;
	for (size_t i=0; i < evolution_ages.size(); i++) 
		Irad_list.push_back(Itot[i]-Iconv[i]);
	Irad = list_to_valarray(Irad_list);

	Mrad_deriv = tabulate_track(exact_track(Mrad_coef, mass),
			evolution_ages, 1);
	Lconv_deriv = tabulate_track(Lconv_track, evolution_ages, 1);
	Lrad_deriv = tabulate_track(Lrad_track, evolution_ages, 1);
	Rrad = tabulate_track(exact_track(Rcore_coef, mass), evolution_ages);
	all_radii = tabulate_track(exact_track(r_coef, mass), evolution_ages);

	std::valarray<double> fake_derivs;
	(*star)->set_angular_momentum_evolution(evolution_ages,
			tabulate_track(Lrad_track, evolution_ages),
			tabulate_track(Lrad_track, evolution_ages, 1), radiative);
	(*star)->set_angular_momentum_evolution(evolution_ages,
			tabulate_track(Lconv_track, evolution_ages),
			tabulate_track(Lconv_track, evolution_ages, 1), convective);
	num_stars_created++;
}

void PlanetData::create_random_planet(Planet** planet)
{
	mass = rand_value(min_planet_mass, max_planet_mass);
	radius = rand_value(min_planet_radius, max_planet_radius);

	std::list<double> semis;
	std::list<double> ages;
	for (double age=0; age < 3; age += 0.1) {
		ages.push_back(age);
		semis.push_back(get_semi(age));
	}
	this->ages.resize(ages.size());
	this->ages = list_to_valarray(ages);
	this->semis.resize(semis.size());
	this->semis = list_to_valarray(semis);

	std::valarray<double> fakeDerivs;
	double curr_semi = get_semi(star->current_age());
	(*planet) = new Planet(star, mass, radius, curr_semi);
	(*planet)->set_semimajor_evolution(this->ages,
			this->semis, fakeDerivs);
}

double PlanetData::get_semi(double age)
{
	return (-age*age*age + 1.5*age*age + 2*age + 3)/50;
}

PlanetData::~PlanetData()
{
	delete sdata;
	if(star) delete star;
}

///Prints an expression of the polynomial that this track represents.
std::ostream &operator<<(std::ostream &os,
			const PolynomialEvolutionTrack &track)
{

	for(size_t p=0; p<track.poly_coef.size(); p++) {
		os << track.poly_coef[p];
		if(p) os << "x";
		if(p>1)	os << "^" << p;
		if(p<track.poly_coef.size()-1) os << " + ";
	}
	os << ", " << track.xmin << " < x < " << track.xmax;
	return os;
}

///Returns an evolutionary track for a quantity that is polynomial in 
///both mass and age, with the given polynomial coefficients.
///The first index should be age and the second one mass.
PolynomialEvolutionTrack *exact_track(
		const std::valarray< std::valarray<double> > &poly_coef, double mass,
		double low_mass_age_scaling, double high_mass_age_scaling,
		double scale_mass)
{
	if(std::isnan(scale_mass)) scale_mass=mass;
	std::valarray<double> age_poly_coeff(poly_coef.size());
	double age_scaling=(mass<=max_low_mass ? low_mass_age_scaling : 
			high_mass_age_scaling);
	for(size_t age_i=0; age_i<poly_coef.size(); age_i++) {
		double mass_pow=1.0;
		age_poly_coeff[age_i]=0.0;
		for(size_t mass_i=0; mass_i<poly_coef[age_i].size(); mass_i++) {
			age_poly_coeff[age_i]+=mass_pow*poly_coef[age_i][mass_i];
			mass_pow*=mass;
		}
		age_poly_coeff[age_i]*=std::pow(scale_mass, age_scaling*age_i);
	}
	return new PolynomialEvolutionTrack(age_poly_coeff, min_age,
			std::numeric_limits<double>::max());
}

///Returns the value of the polynomial with the given coefficients at the
///given mass and age.
double eval_poly(const std::valarray< std::valarray<double> > &poly_coef,
		double mass, double age, double low_mass_age_scaling,
		double high_mass_age_scaling, double scale_mass)
{
	if(std::isnan(scale_mass)) scale_mass=mass;
	double age_scaling=(mass<=max_low_mass ? low_mass_age_scaling : 
			high_mass_age_scaling), result=0.0, age_pow=1.0;
	for(size_t age_i=0; age_i<poly_coef.size(); age_i++) {
		double mass_pow=1.0, mass_result=0.0;
		for(size_t mass_i=0; mass_i<poly_coef[age_i].size(); mass_i++) {
			mass_result+=mass_pow*poly_coef[age_i][mass_i];
			mass_pow*=mass;
		}
		result+=mass_result*age_pow*std::pow(scale_mass, age_scaling*age_i);
		age_pow*=age;
	}
	return result;
}

///Outputs the mass and age polynomial defined by the given polynomial
///coefficients array
std::ostream &operator<<(std::ostream &os, 
		const std::valarray< std::valarray<double> > &poly_coef)
{
	for(size_t age_i=0; age_i<poly_coef.size(); age_i++) {
		for(size_t mass_i=0; mass_i<poly_coef[age_i].size(); mass_i++) {
			os << poly_coef[age_i][mass_i];
			if(age_i) os << "*age";
			if(age_i>1) os << "**" << age_i;
			if(mass_i) os << "*m";
			if(mass_i>1) os << "**" << mass_i;
			if(age_i<poly_coef.size()-1 || mass_i<poly_coef.size()-1)
				os << " + ";
		}
	}
	return os;
}

///Fills the given valarray with a random set of polynomial coefficients.
void rand_poly_coef(std::valarray< std::valarray<double> > &poly_coef,
		double max_mass)
{
	if (max_mass < 0) max_mass = max_stellar_mass;
	poly_coef.resize(3, std::valarray<double>(3));
	double age_fac=1;
	for (unsigned age_i=0; age_i < poly_coef.size(); age_i++) {
		if (age_i > 1) age_fac *= age_i*max_age;
		double mass_fac=1;
		for(unsigned m_i=0; m_i < poly_coef[age_i].size(); m_i++) {
			if(m_i>1) mass_fac*=m_i*max_mass;
			poly_coef[age_i][m_i]=static_cast<double>(rand())/RAND_MAX/
				mass_fac/age_fac;
		}
	}
}

///Returns a random set of polynomial coefficients
std::valarray< std::valarray<double> > rand_poly_coef(double max_mass)
{
	std::valarray< std::valarray<double> > poly_coef(std::valarray<double>(3), 3);
	rand_poly_coef(poly_coef, max_mass);
	return poly_coef;
}

///Returns an array of the values of the track at the given ages.
std::valarray<double> tabulate_track(PolynomialEvolutionTrack *track,
		std::valarray<double> ages, unsigned deriv_order)
{
	std::valarray<double> data(ages.size());
	for(unsigned a=0; a<ages.size(); a++) {
		(*track)(ages[a]);
		data[a]=track->order(deriv_order);
	}
	return data;
}

PolynomialStellarEvolution::PolynomialStellarEvolution(
			const std::valarray<double> &masses,
			const std::valarray<double> &ages,
			const std::valarray< std::valarray<double> > &r_coef,
			const std::valarray< std::valarray<double> > &Iconv_coef,
			const std::valarray< std::valarray<double> > &Itot_coef,
			const std::valarray< std::valarray<double> > &Mrad_coef,
			const std::valarray< std::valarray<double> > &Rcore_coef,
			double low_mass_age_scaling, double high_mass_age_scaling)
{
	std::list< std::valarray<double> > r_data, Iconv_data, Itot_data,
	Mrad_data, Rcore_data;
	PolynomialEvolutionTrack *track;
	for(unsigned i=0; i<masses.size(); i++) {
		track=exact_track(r_coef, masses[i], low_mass_age_scaling,
				high_mass_age_scaling);
		r_data.push_back(tabulate_track(track, ages));
		delete track;
		track=exact_track(Iconv_coef, masses[i],
				low_mass_age_scaling, high_mass_age_scaling);
		Iconv_data.push_back(tabulate_track(track, ages));
		delete track;
		track=exact_track(Itot_coef, masses[i],
				low_mass_age_scaling, high_mass_age_scaling);
		Itot_data.push_back(tabulate_track(track, ages));
		delete track;
		track=exact_track(Mrad_coef, masses[i],
				low_mass_age_scaling, high_mass_age_scaling);
		Mrad_data.push_back(tabulate_track(track, ages));
		delete track;
		track=exact_track(Rcore_coef, masses[i],
				low_mass_age_scaling, high_mass_age_scaling);
		Rcore_data.push_back(tabulate_track(track, ages));
		delete track;
	};
	interpolate_from(masses,
			std::list< std::valarray<double> >(masses.size(), ages),
			r_data, Iconv_data, Itot_data, Mrad_data, Rcore_data, NaN, NaN,
			NaN, std::list< std::valarray<double> >(), max_low_mass,
			low_mass_age_scaling, high_mass_age_scaling);
}

double ExponentialPlusFunc::order(unsigned deriv_order) const
{
	if(deriv_order==0) return (*this)(__deriv_x);
	const FunctionDerivatives *offset_deriv=__offset->deriv(__deriv_x);
	double result=__scale*std::pow(__rate, static_cast<int>(deriv_order))*
		std::exp(__rate*__deriv_x) + offset_deriv->order(deriv_order);
	delete offset_deriv;
	return result;
}

double FuncPlusFunc::order(unsigned deriv_order) const
{
	if(deriv_order==0) return (*this)(__deriv_x);
	const FunctionDerivatives *f1_deriv=__f1->deriv(__deriv_x),
		  *f2_deriv=__f2->deriv(__deriv_x);
	double result=f1_deriv->order(deriv_order)+f2_deriv->order(deriv_order);
	delete f1_deriv;
	delete f2_deriv;
	return result;
}

PiecewiseFunction::PiecewiseFunction(
		const std::list<const OneArgumentDiffFunction *> &pieces,
		double deriv_x) :
	__deriv_x(deriv_x), __range_low(Inf), __range_high(-Inf),__pieces(pieces)
{
	for(std::list<const OneArgumentDiffFunction *>::const_iterator
			piece_i=__pieces.begin(); piece_i!=__pieces.end(); piece_i++) {
		if((*piece_i)->range_low()<__range_low)
			__range_low=(*piece_i)->range_low();
		if((*piece_i)->range_high()>__range_high)
			__range_high=(*piece_i)->range_high();
	}
}

void PiecewiseFunction::add_piece(const OneArgumentDiffFunction *piece)
{
	__pieces.push_back(piece);
	if(piece->range_low()<__range_low)
		__range_low=piece->range_low();
	if(piece->range_high()>__range_high)
		__range_high=piece->range_high();
}

double PiecewiseFunction::operator()(double x) const
{
	double deriv_x=__deriv_x;
	const_cast<PiecewiseFunction*>(this)->__deriv_x=x;
	double result=order(0);
	const_cast<PiecewiseFunction*>(this)->__deriv_x=deriv_x;
	return result;
}

double PiecewiseFunction::order(unsigned deriv_order) const
{
	unsigned index=0;
	if(std::isnan(__deriv_x)) return NaN;
	for(std::list<const OneArgumentDiffFunction *>::const_iterator
			fi=__pieces.begin(); fi!=__pieces.end(); fi++) { 
		if(__deriv_x>=(*fi)->range_low() && __deriv_x<=(*fi)->range_high()) {
			if(deriv_order==0) return (**fi)(__deriv_x);
			const FunctionDerivatives *df=(*fi)->deriv(__deriv_x);
			double result=df->order(deriv_order);
			delete df;
			return result;
		}
		index++;
	}
	std::ostringstream msg;
	msg << "Requested derivative or function value at age=" << __deriv_x
		<< ", outside the range of any piece in PiecewiseFunction::order.";
	throw Error::BadFunctionArguments(msg.str());
}

double FunctionRatio::order(unsigned deriv_order) const
{
	if(deriv_order==0) return (*this)(__deriv_x);
	else {
		const FunctionDerivatives *df1=__f1->deriv(__deriv_x),
			  *df2=__f2->deriv(__deriv_x);
		double result;
		if(deriv_order==1) result=df1->order(1)/df2->order(0) -
			df1->order(0)*df2->order(1)/std::pow(df2->order(0), 2);
		else if(deriv_order==2)
			result=df1->order(2)/df2->order(0) -
				2.0*df1->order(1)*df2->order(1)/std::pow(df2->order(0), 2) -
				df1->order(0)*df2->order(2)/std::pow(df2->order(0), 2) +
				df1->order(0)*std::pow(df2->order(1), 2)/
					std::pow(df2->order(0), 3);
		else throw Error::BadFunctionArguments("Function ratio derivatives "
				"are only implemneted up to and including order 2.");
		delete df1; delete df2;
		return result;
	}
}

double FunctionToPower::order(unsigned deriv_order) const
{
	if(deriv_order==0) return (*this)(__deriv_x);
	else {
		const FunctionDerivatives *df=__f->deriv(__deriv_x);
		double result;
		if(deriv_order==1)
			result=__power*std::pow(df->order(0), __power-1)*df->order(1);
		else if(deriv_order==2)
			result=__power*((__power-1)*std::pow(df->order(0), __power-2)*
					std::pow(df->order(1), 2) +
					std::pow(df->order(0), __power-1)*df->order(2));
		else throw Error::BadFunctionArguments("Function to power "
				"derivatives are only implemneted up to and including order "
				"2.");
		delete df;
		return result;
	}
}

double ScaledFunction::order(unsigned deriv_order) const
{
	if(deriv_order==0) return (*this)(__deriv_x);
	const FunctionDerivatives *f_deriv=__f->deriv(__deriv_x);
	double result=__scale*f_deriv->order(deriv_order);
	delete f_deriv;
	return result;
}

double LogFunction::order(unsigned deriv_order) const
{
	if(deriv_order==0) return (*this)(__deriv_x);
	else {
		const FunctionDerivatives *f_deriv=__f->deriv(__deriv_x);
		double result;
		if(deriv_order==1) result=f_deriv->order(1)/f_deriv->order(0);
		else if(deriv_order==2)
			result=f_deriv->order(2)/f_deriv->order(0) -
				std::pow(f_deriv->order(1)/f_deriv->order(0), 2);
		else throw Error::BadFunctionArguments("Log(Function) derivatives "
				"are only implemneted up to and including order 2.");
		delete f_deriv;
		return result;
	}
}

double rand_value(double min, double max)
{
	return (max-min)*rand()/RAND_MAX+min;
}

int rand_value(int min, int max)
{
	return rand()%(max-min+1)+min;
}

///Given n, m and (n)C(m) return (n)C(m+1)
unsigned next_binom_coef(unsigned n, unsigned m, unsigned nCm)
{
	assert(n>=m+1);
	return (nCm*(n-m))/(m+1);
}

///Returns new polynomial coefficienst such that
///output polynomial(mass, age+age_offset)=input polynomial(mass, age)
std::valarray< std::valarray<double> > offset_age(
		const std::valarray< std::valarray<double> > &poly_coef,
		double age_offset)
{
	std::valarray< std::valarray<double> > result(
			std::valarray<double>(poly_coef[0].size()), poly_coef.size());
	for(size_t i=0; i<poly_coef.size(); i++) 
		for(size_t j=0; j<poly_coef[0].size(); j++) {
			double c_ij=poly_coef[i][j], offset_k=1.0;
			size_t binom_coef=1;
			for(size_t k=0; k<=i; k++) {
				result[i-k][j]+=c_ij*offset_k*binom_coef;
				offset_k*=age_offset;
				if(k<i) binom_coef=next_binom_coef(i, k, binom_coef);
			}
		}
	return result;
}

///Solves the given equation and its derivative using
double solve(double guess_x, double abs_precision, double rel_precision,
		double (*f)(double x, void *params),
		double (*df) (double x, void *params),
		void (*fdf) (double x, void *params, double *f, double *df),
		void *params)
{
	int status;
	int iter = 0;
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;
	double x0, x = guess_x;
	gsl_function_fdf FDF;

	FDF.f = f;
	FDF.df = df;
	FDF.fdf = fdf;
	FDF.params = params;

	T = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc (T);
	gsl_root_fdfsolver_set (s, &FDF, x);

	do {
		iter++;
		status = gsl_root_fdfsolver_iterate (s);
		x0 = x;
		x = gsl_root_fdfsolver_root (s);
		status = gsl_root_test_delta (x, x0, abs_precision, rel_precision);
	} while (status == GSL_CONTINUE);

	gsl_root_fdfsolver_free (s);
	return x;
}
