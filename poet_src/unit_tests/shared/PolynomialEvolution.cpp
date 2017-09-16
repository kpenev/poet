#include "PolynomialEvolution.h"

namespace StellarEvolution {

    double PolynomialEvolutionQuantity::order(unsigned deriv_order) const
    {
        if(__deriv_x < __xmin || __deriv_x > __xmax) {
            std::ostringstream msg;
            msg << ("Attempt to evaluate a PolynomialEvolutionQuantity "
                    "defined over ")
                << __xmin << " < x < " << __xmax << " at x=" << __deriv_x
                << ", polynomial coefficients:" << __poly_coef;
            std::cout << msg.str() << std::endl;
//            *static_cast<int*>(NULL) = 0;
            throw Core::Error::Runtime(msg.str());
        }
        double xn = 1, result = 0;
        for(size_t n = deriv_order; n < __poly_coef.size(); ++n) {
            double coef_mod = 1;
            for(size_t i = 0; i < deriv_order; ++i) coef_mod *= (n-i);
            result += coef_mod * __poly_coef[n] * xn;
            xn *= __deriv_x;
        }
        return result;
    }

    MockStellarEvolution::MockStellarEvolution(
        double core_formation_age,
        const std::valarray< std::valarray<double> > &R,
        const std::valarray< std::valarray<double> > &Iconv,
        const std::valarray< std::valarray<double> > &Irad,
        const std::valarray< std::valarray<double> > &Rcore,
        const std::valarray< std::valarray<double> > &Mcore,
        const std::valarray< std::valarray<double> > &Lum
    ) :
        __R(R),
        __Iconv(Iconv),
        __Irad(Irad),
        __Rcore(Rcore),
        __Mcore(Mcore),
        __Lum(Lum),
        __core_formation_age(core_formation_age)
    {
        if(std::isnan(__core_formation_age))
            __core_formation_age = rand_value(MIN_AGE, 2.1);
        if(__R.size() == 0) rand_poly_coef(__R);
        if(__Iconv.size()==0) rand_poly_coef(__Iconv);
        if(__Irad.size()==0) rand_poly_coef(__Irad);
        if(__Rcore.size()==0) rand_poly_coef(__Rcore);
        if(__Mcore.size()==0) rand_poly_coef(__Mcore);
        if(__Lum.size()==0) rand_poly_coef(__Lum);
        if(core_formation_age > 0) {
            __Irad[0] = 0;
            __Rcore[0] = 0;
            __Mcore[0] = 0;
            __Irad = offset_age(Irad, core_formation_age);
            __Rcore = offset_age(Rcore, core_formation_age);
            __Mcore = offset_age(Mcore, core_formation_age);
        }
        __Menv.resize(
            __Mcore.size(),
            std::valarray<double>(0.0,
                                  std::max(__Mcore[0].size(), size_t(2)))
        );
        __Itot.resize(__Iconv.size(),
                      std::valarray<double>(__Iconv[0].size()));
        __Itot = __Iconv + __Irad;
        for(size_t mass_i = 0; mass_i < __Mcore.size(); ++mass_i)
            for(size_t age_i = 0; age_i < __Mcore[mass_i].size(); ++age_i)
                __Menv[mass_i][age_i] = -__Mcore[mass_i][age_i];
        __Menv[0][1] += 1;
    }

    EvolvingStellarQuantity *MockStellarEvolution::operator()(
        QuantityID quantity,
        double mass,
        double
    ) const
    {
        switch(quantity) {
            case RADIUS: return exact_track(__R, mass);
            case ICONV: return exact_track(__Iconv, mass);
            case LUM: return exact_track(__Lum, mass);
            case IRAD: return exact_track(__Irad, mass);
            case MRAD: return exact_track(__Mcore, mass);
            case RRAD: return exact_track(__Rcore, mass);
            default: assert(false);
        };
    }
    
}//End StellarEvolution namespace.

#if 0
    StarData::StarData() :
        num_stars_created(0),
        Lrad_track(NULL),
        Lconv_track(NULL),
        evolution_masses(100),
        evolution_ages(100),
        r_coef(rand_poly_coef()),
        Iconv_coef(rand_poly_coef()),
        Irad_coef(rand_poly_coef()),
        Mrad_coef(rand_poly_coef()),
        Rcore_coef(rand_poly_coef())
    {
        for (size_t age_ind=0; age_ind < evolution_ages.size(); age_ind++) {
            evolution_ages[age_ind] = (
                MIN_AGE
                +
                age_ind * (max_age - MIN_AGE) / evolution_ages.size()
            );
        }

        for (
            size_t mass_ind = 0;
            mass_ind < evolution_masses.size();
            ++mass_ind
        ) {
            evolution_masses[mass_ind] = (
                min_stellar_mass
                +
                mass_ind
                *
                (max_stellar_mass - min_stellar_mass)
                /
                evolution_masses.size()
            );
        }
    }


    StarData::~StarData()
    {
        if(Lrad_track) delete Lrad_track;
        if(Lconv_track) delete Lconv_track;
    }

    void StarData::create_random_star(Star::InterpolatedEvolutionStar **star)
    {
        mass = rand_value(min_stellar_mass, max_stellar_mass);
        feh = rand_value(-1.0, 0.5);
        age = rand_value(MIN_AGE, 2.1);
        radius = (rand_value(0.9, 1.1)) * eval_poly(r_coef, mass, age);
        conv_spin = rand_value(0.0, 1.0);
        rad_spin = rand_value(0.0, 1.0);
        tidal_Q = std::pow(10.0, rand_value(5.0,10.0));
        wind_strength = rand_value(0.0, 1.0);
        wind_sat_freq = rand_value(0.0, 1.0);
        disk_lock_w = rand_value(0.0, wind_sat_freq);
        disk_lock_time = evolution_ages[rand_value(0, evolution_ages.size()-1)];
        coupling_timescale = rand_value(0.0, 1.0);
        Q_trans_width = 0;

        MockStellarEvolution evol(Core::NaN,
                                  r_coef,
                                  Iconv_coef,
                                  Irad_coef,
                                  Rcore_coef,
                                  Mcore_coef);
        (*star) = new Star(mass,
                           feh,
                           wind_strength,
                           wind_sat_freq,
                           coupling_timescale,
                           evol);

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
	double curr_semi = get_semi(star->age());
	(*planet) = new Planet(*star, mass, radius, curr_semi);
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
#endif

namespace StellarEvolution {

    PolynomialEvolutionQuantity *exact_track(
        const std::valarray< std::valarray<double> > &poly_coef,
        double mass,
        double low_mass_age_scaling,
        double high_mass_age_scaling,
        double scale_mass
    )
    {
        if(std::isnan(scale_mass)) scale_mass=mass;
        std::valarray<double> age_poly_coeff(poly_coef.size());
        double age_scaling=(mass <= MAX_LOW_MASS ? low_mass_age_scaling
                                                 : high_mass_age_scaling);
        for(size_t age_i = 0; age_i < poly_coef.size(); ++age_i) {
            double mass_pow=1.0;
            age_poly_coeff[age_i]=0.0;
            for(size_t mass_i=0; mass_i<poly_coef[age_i].size(); mass_i++) {
                age_poly_coeff[age_i]+=mass_pow*poly_coef[age_i][mass_i];
                mass_pow*=mass;
            }
            age_poly_coeff[age_i]*=std::pow(scale_mass, age_scaling*age_i);
        }
        return new PolynomialEvolutionQuantity(
            age_poly_coeff,
            MIN_AGE,
            std::numeric_limits<double>::max()
        );
    }

}//End StellarEvolution namespace.

#if 0
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

std::valarray<double> tabulate_track(PolynomialEvolutionQuantity *track,
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
	std::list< std::valarray<double> > r_data, Iconv_data, Irad_data,
	Mrad_data, Rcore_data;
	PolynomialEvolutionQuantity *track;
	for(unsigned i=0; i<masses.size(); i++) {
		track=exact_track(r_coef, masses[i], low_mass_age_scaling,
				high_mass_age_scaling);
		r_data.push_back(tabulate_track(track, ages));
		delete track;
		track=exact_track(Iconv_coef, masses[i],
				low_mass_age_scaling, high_mass_age_scaling);
		Iconv_data.push_back(tabulate_track(track, ages));
		delete track;
		track=exact_track(Itot_coef-Iconv_coef, masses[i],
				low_mass_age_scaling, high_mass_age_scaling);
		Irad_data.push_back(tabulate_track(track, ages));
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
			r_data, Iconv_data, Irad_data, Mrad_data, Rcore_data, NaN, NaN,
			NaN, NaN, NaN, 0, 0, 0, 0, 0, 
			std::list< std::valarray<double> >(), 0,
			max_low_mass, low_mass_age_scaling, high_mass_age_scaling);
}
#endif

double ExponentialPlusFunc::order(unsigned deriv_order) const
{
	if(deriv_order == 0) return (*this)(__deriv_x);
	const Core::FunctionDerivatives *offset_deriv = __offset->deriv(__deriv_x);
	double result = (__scale
                     *
                     std::pow(__rate, static_cast<int>(deriv_order))
                     *
                     std::exp(__rate * __deriv_x)
                     +
                     offset_deriv->order(deriv_order));
	delete offset_deriv;
	return result;
}

#if 0
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
#endif
