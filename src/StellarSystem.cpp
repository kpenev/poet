/**\file
 *
 * \brief The definition of some of the methods of the StellarSystem class.
 *
 * \ingroup StellarSystem_group
 */

#include "StellarSystem.h"
#include "Error.h"
#include "OrbitSolver.h"
#include <limits>
#include <iostream>

#include <gsl/gsl_matrix.h>

/*
double energy(void* xp) {
	StellarSystem* system = static_cast< StellarSystem* >(xp);
	return system->annealing_energy();
}

void step(const gsl_rng *r, void *xp, double step_size) {
	StellarSystem* system = static_cast< StellarSystem* >(xp);
	system->annealing_step(r, step_size);
}

double metric(void *xp, void *yp) {
	StellarSystem* first = static_cast< StellarSystem* >(xp);
	StellarSystem* second = static_cast< StellarSystem* >(yp);
	return first->annealing_metric(second);
}

void print(void* xp) {
	StellarSystem* system = static_cast< StellarSystem* >(xp);
	system->annealing_print();
}

bool stop_evolution(double age, const double* y, void *stellar_system) {
	if (std::isnan(y[0]) || std::isnan(y[1]) || std::isnan(y[2]))
		return true;
	StellarSystem* system =
			static_cast< StellarSystem* > (stellar_system);
	//if (age >= system->get_star().get_lifetime()) return true;
	double min_semi = system->get_planet().minimum_semimajor(age);
	double curr_semi = std::pow(y[0], 1.0/6.5)*
			AstroConst::solar_radius/AstroConst::AU;
	//using namespace std;
	//cout<<"stop evol "<<curr_semi<<" "<<min_semi<<endl;
	if (std::isnan(curr_semi)) return true;
	return curr_semi <= min_semi;
}*/

StellarSystem::StellarSystem(Star *system_star, Planet *system_planet,
		double age, const std::string &system_name) :
	name(system_name), age(age), star(system_star), planet(system_planet)
{
}

void StellarSystem::transform_into_earth() {
	planet->transform_into_earth();
}

//unused for now
/*
///Make a stellar system from the given star and planet having the
///given age and name.
StellarSystem::StellarSystem(
		Star *system_star,
		Planet *system_planet,
		OrbitSolver solve_ode) :
	name(system_star->get_sysname()), age(system_star->current_age()),
	star(system_star), planet(system_planet)
{
	using namespace std;
	using namespace AstroConst;
	double present_orbit[3]={
		std::pow(planet->get_current_semimajor(), 6.5), 
		star->get_current_angular_momentum(convective),
		star->get_current_angular_momentum(radiative)};
	present_orbit[0] *= std::pow(AstroConst::AU/AstroConst::solar_radius,6.5);
	solve_ode(stellar_system_diff_eq, stellar_system_jacobian,
			age, present_orbit, stop_evolution, (void*)this);
	const std::list<double> 
		*evol_ages=solve_ode.get_tabulated_indep_var(),
		*evol_a6p5=solve_ode.get_tabulated_dep_var(0),
		*evol_Lconv=solve_ode.get_tabulated_dep_var(1),
		*evol_Lrad=solve_ode.get_tabulated_dep_var(2),
		*evol_da6p5_dt=solve_ode.get_tabulated_dep_var_deriv(0),
		*evol_dLconv_dt=solve_ode.get_tabulated_dep_var_deriv(1),
		*evol_dLrad_dt=solve_ode.get_tabulated_dep_var_deriv(2);
	std::valarray<double> age_array=list_to_valarray(*evol_ages), 
						  a_array(evol_a6p5->size()),
						  da_dt_array(evol_a6p5->size());
	std::list<double>::const_iterator a6p5=evol_a6p5->begin(), 
									  da6p5_dt=evol_da6p5_dt->begin();
	if(evol_ages->size()!=evol_a6p5->size() ||
					evol_ages->size()!=evol_da6p5_dt->size())
		throw Error::Runtime("The evolution of a^6.5 or d(a^6.5)/dt does not"
			" contain the same number of points as the tabulated ages. in"
			"StellarSystem constructor.");
	for(size_t i=0; i<evol_ages->size(); i++) {
		a_array[i]=std::pow(*a6p5, 1.0/6.5)*Rsun_AU;
		da_dt_array[i]=*da6p5_dt/(6.5*std::pow(*a6p5, 11.0/13.0))*Rsun_AU;
		da6p5_dt++;
		a6p5++;

	}
	std::valarray<double> Lconv_arr = list_to_valarray(*evol_Lconv),
			Lrad_arr = list_to_valarray(*evol_Lrad);
	for (size_t i=0; i < evol_ages->size(); i++) {
		std::cout<<age_array[i]<<" "<<a_array[i]<<" "<<Lconv_arr[i]<<" "<<
				Lrad_arr[i]<<std::endl;
		double Lorb = planet->orbital_angmom(a_array[i]);
		cout << "angmom "<<Lconv_arr[i]+Lrad_arr[i]+Lorb<<endl;

	}
	try{
	planet->set_semimajor_evolution(age_array, a_array, da_dt_array);
	star->set_angular_momentum_evolution(age_array,
			list_to_valarray(*evol_Lconv), list_to_valarray(*evol_dLconv_dt),
			convective);
	star->set_angular_momentum_evolution(age_array,
			list_to_valarray(*evol_Lrad), list_to_valarray(*evol_dLrad_dt),
			radiative);
	}
	catch(alglib::ap_error &e) {
		cout<<e.msg<<endl;
	}
}*/

/*
double StellarSystem::annealing_energy() {
	double AU_Rsun = AstroConst::AU/AstroConst::solar_radius;
	OrbitSolver solver(MIN_AGE, age, PRECISION);
	star->set_disk_lock_frequency(w0);
	solver(*this, Inf, AU_Rsun*a0);
	double frac_error = solver.last_error(age, planet->get_current_semimajor(),
			star->get_current_angular_momentum(convective));
	if (std::isnan(frac_error)) return INVALID_ENERGY;
	return std::log10(frac_error);
}

void StellarSystem::annealing_step(const gsl_rng *r, double step_size) {
	double mid_a0 = (min_a0 + max_a0)/2;
	double mid_w0 = (min_w0 + max_w0)/2;
	//double step_size = max_step;
	step_size *= gsl_rng_uniform(r);
	double new_a0 = a0;
	double new_w0 = w0;
	while (new_a0 == a0 && new_w0 == w0) {
		double angle = gsl_rng_uniform(r)*2*M_PI;
		new_a0 += mid_a0*step_size*std::cos(angle);
		new_w0 += mid_w0*step_size*std::sin(angle);
		new_a0 = std::min(max_a0, new_a0);
		new_a0 = std::max(min_a0, new_a0);
		new_w0 = std::min(max_w0, new_w0);
		new_w0 = std::max(min_w0, new_w0);
	}
	a0 = new_a0;
	w0 = new_w0;
}

double StellarSystem::annealing_metric(StellarSystem* that) {
	double mid_a0 = (min_a0 + max_a0)/2;
	double mid_w0 = (min_w0 + max_w0)/2;
	double dist = std::pow((this->a0 - that->a0)/mid_a0, 2) +
			std::pow((this->w0 - that->w0)/mid_w0, 2);
	return std::sqrt(dist);
}

void StellarSystem::annealing_print() {
	std::cout<<" a0 w0 "<<a0<<" "<<w0<<" ";
}

bool StellarSystem::anneal_solve_IC(bool verbose, double min_a0,
		double max_a0, double min_w0, double max_w0, double max_err) {
	this->min_a0 = min_a0;
	this->max_a0 = max_a0;

	this->min_w0 = min_w0;
	this->max_w0 = max_w0;
	a0 = (min_a0 + max_a0)/2;
	w0 = (min_w0 + max_w0)/2;

	gsl_siman_params_t params;
	params.iters_fixed_T = 2;
	params.k = 0.2;
	params.mu_t = 1.01;
	params.n_tries = 4;
	params.t_initial = 0.1;
	params.t_min = 1e-3;
	params.step_size = 2e-2;
	gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (r, time(NULL));
	gsl_siman_print_t print_func;
	if (verbose) print_func = print;
	else print_func = NULL;
	gsl_siman_solve (r, (void*) this, energy, step, metric, print_func,
			NULL, NULL, NULL, sizeof(StellarSystem), params);
	gsl_rng_free(r);
	double last_energy = annealing_energy();
	std::cout<<name<<" "<<a0<<" "<<w0<<" "<<last_energy<<std::endl;
	if (last_energy == INVALID_ENERGY || last_energy > std::log10(max_err))
		return false;
	return true;
}*/

int StellarSystem::orbit_differential_equation(
		double age,
		const double* orbital_parameters,
		double* orbital_derivatives, int assume_sign,
		WindSaturationState assume_wind_saturation) const
{
	double wconv=star->spin_frequency(age, convective,
			orbital_parameters[1]);
	if(orbital_parameters[0]<0) return GSL_ERANGE;
	double semimajor=std::pow(std::abs(orbital_parameters[0]), 1.0/6.5)*
		Rsun_AU;
	if(assume_sign==0) 
		orbital_derivatives[0]=planet->tidal_decay(age,semimajor,wconv,true);
	else if(assume_sign>0)
		orbital_derivatives[0]=planet->tidal_decay(age, semimajor, Inf,true);
	else
		orbital_derivatives[0]=planet->tidal_decay(age, semimajor, -Inf,
				true);
	orbital_derivatives[0] /= std::pow(Rsun_AU, 6.5);
	double semi_deriv = orbital_derivatives[0]/
			(6.5*std::pow(orbital_parameters[0], 11.0/13.0))*Rsun_AU;
	no_planet_differential_equation(age, orbital_parameters+1,
			orbital_derivatives+1, assume_wind_saturation);
	orbital_derivatives[1]-=planet->orbital_angmom_deriv(semimajor,
			semi_deriv);
	return GSL_SUCCESS;
}

int StellarSystem::no_planet_differential_equation(double age,
		const double *orbital_parameters, double *orbital_derivatives,
		WindSaturationState assume_wind_saturation) const
{
	double coupling_torque=star->differential_rotation_torque_angmom(age,
			orbital_parameters[0], orbital_parameters[1]),
		   wconv=star->spin_frequency(age, convective,orbital_parameters[0]);
	orbital_derivatives[0]=coupling_torque -
		star->wind_torque(age, wconv, assume_wind_saturation);
	orbital_derivatives[1]=-coupling_torque;
	return GSL_SUCCESS;
}

int StellarSystem::no_planet_jacobian(double age,
		const double *orbital_parameters, double *param_derivs,
		double *age_derivs,
		WindSaturationState assume_wind_saturation) const
{
	double Lconv=orbital_parameters[0],
		   Lrad=orbital_parameters[1],
		   wconv=star->spin_frequency(age, convective, Lconv),
		   dwconv_dLconv=star->spin_frequency_angmom_deriv(age, convective,
				Lconv),
		   dcoupling_torque_dLconv=star->differential_rotation_torque_deriv(
				   age, Lconv, Lrad, convective),
		   dcoupling_torque_dLrad=star->differential_rotation_torque_deriv(
				   age, Lconv, Lrad, radiative),
		   dcoupling_torque_dage=star->differential_rotation_torque_deriv(
				   age, Lconv, Lrad);
	param_derivs[2]=-dcoupling_torque_dLconv;
	param_derivs[3]=-dcoupling_torque_dLrad;
	param_derivs[0]=-param_derivs[2]
		-star->wind_torque_freq_deriv(age, wconv, assume_wind_saturation)*
		dwconv_dLconv;
	param_derivs[1]=-param_derivs[3];
	age_derivs[0]=dcoupling_torque_dage -
		star->wind_torque_age_deriv(age, Lconv, true,
				assume_wind_saturation);
	age_derivs[1]=-dcoupling_torque_dage;
	return GSL_SUCCESS;
}

///The differential equation governing the rotation of the star's
///radiative zone with the convective zone locked to a disk: takes age
///and pointer to the radiative antgular momentum (in units of
///Msun*Rsun^2/day) dependent variables on input and updates the
///derivative of Lrad on output.
int StellarSystem::locked_conv_differential_equation(double age,
		const double *orbital_parameters, double *orbital_derivatives)
{
	orbital_derivatives[0]=-star->differential_rotation_torque_angmom(age,
			star->get_disk_lock_frequency()*
			star->moment_of_inertia(age, convective), orbital_parameters[0]);
	return GSL_SUCCESS;
}

///The jacobian of the differential equation governing the rotation of
///the star's radiative zone with the convective zone locked to a disk:
///takes age and pointer to the radiative antgular momentum (in units of
///Msun*Rsun^2/day) dependent variables on input and updates the jacobian
///and partial time derivative of the Lrad equation on output.
int StellarSystem::locked_conv_jacobian(double age, const double* orbital_parameters,
		double* param_derivs, double *age_derivs) const
{
	double wconv=star->get_disk_lock_frequency(),
		   Lconv=wconv*star->moment_of_inertia(age, convective),
		   Lrad=orbital_parameters[0];
	param_derivs[0]=-star->differential_rotation_torque_deriv(
			age, Lconv, Lrad, radiative),
	age_derivs[0] = -star->differential_rotation_torque_deriv(
				   age, Lconv, Lrad, convective)*wconv*
		star->moment_of_inertia_deriv(age, convective) -
		star->differential_rotation_torque_deriv(age, Lconv, Lrad);
	return GSL_SUCCESS;
}

int StellarSystem::orbit_jacobian(double age,
		const double* orbital_parameters, double* param_derivs, 
		double* age_derivs, int assume_sign,
		WindSaturationState assume_wind_saturation) const
{
	double semimajor=std::pow(std::abs(orbital_parameters[0]), 1.0/6.5)*
		Rsun_AU,
		   Lconv=orbital_parameters[1],
		   Lrad=orbital_parameters[2],
		   wconv=star->spin_frequency(age, convective, Lconv),
		   dwconv_dLconv=star->spin_frequency_angmom_deriv(age, convective,
				Lconv),
		   dcoupling_torque_dLconv=star->differential_rotation_torque_deriv(
				   age, Lconv, Lrad, convective),
		   dcoupling_torque_dLrad=star->differential_rotation_torque_deriv(
				   age, Lconv, Lrad, radiative),
		   dcoupling_torque_dage=star->differential_rotation_torque_deriv(
				   age, Lconv, Lrad),
		   wconv_tidal=(assume_sign==0 ? wconv :
				   (assume_sign>0 ? Inf : -Inf)),
		   semi_deriv = planet->tidal_decay(age, semimajor, wconv_tidal,
				   false);
	param_derivs[0]= planet->tidal_decay_semimajor_deriv(age, semimajor,
			wconv_tidal, true);
	param_derivs[2]=0.0;
	param_derivs[1]=planet->tidal_decay_star_spin_deriv(age, semimajor,
			wconv_tidal, true)*-dwconv_dLconv/std::pow(Rsun_AU, 6.5);

	param_derivs[6]=0;
	param_derivs[7]=-dcoupling_torque_dLconv;
	param_derivs[8]=-dcoupling_torque_dLrad;

	param_derivs[3]=-planet->orbital_angmom_deriv_semimajor_deriv(semimajor,
			semi_deriv, wconv_tidal)/(6.5*std::pow(semimajor, 5.5))*
			std::pow(Rsun_AU, 6.5);
	param_derivs[4]=-param_derivs[7]
		-star->wind_torque_freq_deriv(age, wconv, assume_wind_saturation)*
		dwconv_dLconv +
		planet->tidal_torque_star_spin_deriv(age, semimajor, semi_deriv,
				wconv_tidal)*
		-dwconv_dLconv;
	param_derivs[5]=-param_derivs[8];

	age_derivs[0]=planet->tidal_decay_age_deriv(age, semimajor, wconv_tidal,
			true)/pow(Rsun_AU, 6.5);
	age_derivs[2]= -dcoupling_torque_dage;
	age_derivs[1]=-age_derivs[2] -
		star->wind_torque_age_deriv(age, Lconv, true, assume_wind_saturation)
			-planet->orbit_angmom_deriv_age_deriv(age, semimajor, semi_deriv,
					wconv_tidal);
	return GSL_SUCCESS;
}

/*
int StellarSystem::orbit_diffeq_prescribed_rot(
		double age,
		const double* orbital_parameters,
		double* orbital_derivatives) const
{
	using namespace std;
	double wconv=star->spin_frequency(age, convective,
			orbital_parameters[1]),
		   coupling_torque=star->differential_rotation_torque_angmom(age,
				   orbital_parameters[1], orbital_parameters[2]),
		   semimajor=std::pow(orbital_parameters[0], 1.0/6.5)*Rsun_AU;
	orbital_derivatives[0]=planet->tidal_decay(age, semimajor, wconv, true);
	orbital_derivatives[0] /= pow(Rsun_AU, 6.5);
	orbital_derivatives[1] = star->angular_momentum_deriv(age, convective, 1);
	orbital_derivatives[2] = star->angular_momentum_deriv(age, radiative, 1);
	return GSL_SUCCESS;
}

int StellarSystem::orbit_jacobian_prescribed_rot(
		double age,
		const double* orbital_parameters,
		double* param_derivs,
		double* age_derivs) const
{
	double semimajor=std::pow(orbital_parameters[0], 1.0/6.5)*Rsun_AU,
		   Lconv=orbital_parameters[1],
		   Lrad=orbital_parameters[2],
		   wconv=star->spin_frequency(age, convective, Lconv),
		   dwconv_dLconv=star->spin_frequency_angmom_deriv(age, convective,
				Lconv);
	double semi_deriv = planet->tidal_decay(age, semimajor, wconv, false);
	param_derivs[0]= planet->tidal_decay_semimajor_deriv(age, semimajor, wconv, true);
	param_derivs[2]=0.0;
	param_derivs[1]=planet->tidal_decay_star_spin_deriv(age, semimajor,
			wconv, true)*-dwconv_dLconv/std::pow(Rsun_AU, 6.5);

	param_derivs[6]=0;
	param_derivs[7]=0;
	param_derivs[8]=0;

	param_derivs[3]=0;
	param_derivs[4]=0;
	param_derivs[5]=0;

	age_derivs[0]=planet->tidal_decay_age_deriv(age, semimajor, wconv, true)/
			pow(Rsun_AU, 6.5);
	age_derivs[1] = star->angular_momentum_deriv(age, convective, 2);
	age_derivs[2] = star->angular_momentum_deriv(age, radiative, 2);
	return GSL_SUCCESS;
}*/

int StellarSystem::locked_orbit_differential_equation(double age,
		const double* orbital_parameters, double* orbital_derivatives,
		WindSaturationState assume_wind_saturation) const
{
	if(orbital_parameters[0]<0) return GSL_ERANGE;
	double semimajor_rsun=orbital_parameters[0],
		   semimajor_meters=semimajor_rsun*AstroConst::solar_radius,
		   semimajor_au=semimajor_rsun*Rsun_AU,
		   wconv=planet->orbital_angular_velocity_semimajor(semimajor_au),
		   Iconv=star->moment_of_inertia(age, convective),
		   dIconv_dt=star->moment_of_inertia_deriv(age, convective),
		   Lconv=wconv*Iconv,
		   coupling_torque=star->differential_rotation_torque_angmom(age,
				   Lconv, orbital_parameters[1]),
		   wind_torque=star->wind_torque(age, wconv, assume_wind_saturation),
		   mstar=star->get_mass()*AstroConst::solar_mass,
		   mplanet=planet->get_mass()*AstroConst::jupiter_mass,
		   mtotal=mstar+mplanet,
		   torque_coef=std::sqrt(semimajor_meters/mtotal/AstroConst::G)*
			   AstroConst::solar_radius/AstroConst::day;
	orbital_derivatives[0]=2.0*(torque_coef*(coupling_torque-wind_torque)
		- dIconv_dt/semimajor_rsun)/
		(mstar*mplanet/(mtotal*AstroConst::solar_mass) -
		 3.0*Iconv/std::pow(semimajor_rsun,2));
	orbital_derivatives[1]=-coupling_torque;
	return GSL_SUCCESS;
}

int StellarSystem::locked_orbit_jacobian(double age,
		const double* orbital_parameters, double* param_derivs,
		double* age_derivs,
		WindSaturationState assume_wind_saturation) const
{
	double semimajor_rsun=orbital_parameters[0],
		   semimajor_meters=semimajor_rsun*AstroConst::solar_radius,
		   semimajor_au=semimajor_meters/AstroConst::AU,
		   wconv=planet->orbital_angular_velocity_semimajor(semimajor_au),
		   Iconv=star->moment_of_inertia(age, convective),
		   dIconv_dt=star->moment_of_inertia_deriv(age, convective),
		   d2Iconv_dt2=star->moment_of_inertia_deriv(age, convective, 2),
		   Lconv=wconv*Iconv,
		   Lrad=orbital_parameters[1],
		   dwconv_da=planet->orbital_angular_velocity_semimajor_deriv(
				   semimajor_au)*Rsun_AU,
		   dLconv_da=dwconv_da*Iconv,
		   dLconv_dt=wconv*dIconv_dt,
		   coupling_torque=star->differential_rotation_torque_angmom(age,
				   Lconv, Lrad),
		   dcoupling_torque_dLconv=star->differential_rotation_torque_deriv(
				   age, Lconv, Lrad, convective),
		   dcoupling_torque_dLrad=star->differential_rotation_torque_deriv(
				   age, Lconv, Lrad, radiative),
		   dcoupling_torque_da=dcoupling_torque_dLconv*dLconv_da,
		   dcoupling_torque_dt=dcoupling_torque_dLconv*dLconv_dt +
			   star->differential_rotation_torque_deriv(age, Lconv, Lrad),
		   wind_torque=star->wind_torque(age, wconv, assume_wind_saturation),
		   dwind_torque_da=star->wind_torque_freq_deriv(age, wconv,
				   assume_wind_saturation)*dwconv_da,
		   dwind_torque_dt=star->wind_torque_age_deriv(age, Lconv, false,
				   assume_wind_saturation),
		   total_torque=(coupling_torque-wind_torque),
		   mstar=star->get_mass()*AstroConst::solar_mass,
		   mplanet=planet->get_mass()*AstroConst::jupiter_mass,
		   mtotal=mstar+mplanet,
		   torque_coef=std::sqrt(semimajor_meters/mtotal/AstroConst::G)/
			   AstroConst::day*AstroConst::solar_radius,
		   dtorque_coef_da=0.5*std::pow(AstroConst::solar_radius,2)/
			   AstroConst::day/
			   std::sqrt(semimajor_meters*mtotal*AstroConst::G),
		   denominator=(mstar*mplanet/(mtotal*AstroConst::solar_mass) -
				   3.0*Iconv/std::pow(semimajor_rsun,2));

	param_derivs[0]=2.0*(dtorque_coef_da*total_torque +
			torque_coef*(dcoupling_torque_da-dwind_torque_da) + 
			dIconv_dt/std::pow(semimajor_rsun, 2))/denominator - 
		12.0*(torque_coef*total_torque - dIconv_dt/semimajor_rsun)/
		std::pow(denominator, 2)*Iconv/std::pow(semimajor_rsun, 3);

	param_derivs[1]=2.0*(torque_coef*dcoupling_torque_dLrad)/denominator;

	param_derivs[2]=-dcoupling_torque_da;
	param_derivs[3]=-dcoupling_torque_dLrad;

	age_derivs[0]=2.0*(torque_coef*(dcoupling_torque_dt-dwind_torque_dt) - 
			d2Iconv_dt2/semimajor_rsun)/denominator + 
		6.0*(torque_coef*total_torque - dIconv_dt/semimajor_rsun)/
		std::pow(denominator, 2)*dIconv_dt/std::pow(semimajor_rsun, 2);

	age_derivs[1]=-dcoupling_torque_dt;

	return GSL_SUCCESS;
}

/*
void StellarSystem::solve_init(double start_age, double curr_age,
		double curr_spin) {
	using namespace std;
	const int INIT_SIZE = 10;
	const double DEFAULT_INIT_SPIN = 3;
	valarray<double> ages(INIT_SIZE);
	valarray<double> Lc(INIT_SIZE), Lr(INIT_SIZE);
	valarray<double> Lc_deriv(INIT_SIZE), Lr_deriv(INIT_SIZE);
	for (int i=0; i < INIT_SIZE; i++) {
		double age = i*curr_age/INIT_SIZE;
		ages[i] = age;
		Lc[i] = Lr[i] = 0;
		Lc_deriv[i] = Lr_deriv[i] = 0;
	}
	star->set_angular_momentum_evolution(ages, Lc, Lc_deriv, convective);
	star->set_angular_momentum_evolution(ages, Lr, Lr_deriv, radiative);
	double init_Ic = star->moment_of_inertia(start_age, convective);
	double init_Ir = star->moment_of_inertia(start_age, radiative);
	double curr_Ic = star->moment_of_inertia(curr_age, convective);

	OrbitSolver backSolver(curr_age, start_age, PRECISION);
	const double AU_Rsun = AstroConst::AU/AstroConst::solar_radius;
	while (1) {
		//OrbitSolver solver(start_age, curr_age, precision);
		double orbit[3]={
			std::pow(AU_Rsun*planet->get_current_semimajor(), 6.5),
			star->get_angular_momentum(curr_age, convective),
			star->get_angular_momentum(curr_age, radiative)};
		backSolver(system_diffeq_prescribed, system_jacobian_prescribed,
				curr_age, orbit, stop_evolution, (void*)this);
		double init_semi = backSolver.get_tabulated_dep_var(0)->back();
		cout<<"init semi "<<pow(init_semi,1.0/6.5)/AU_Rsun<<endl;
		double init_Lc = DEFAULT_INIT_SPIN*init_Ic;
		double init_Lr = DEFAULT_INIT_SPIN*init_Ir;
		//orbit[0] =
		OrbitSolver forwardSolver(start_age, curr_age, PRECISION);
		cout<<"outer loop"<<endl;
		double min_diff = std::numeric_limits<double>::infinity();
		double best_next_Lc = -1;
		while (1 {
			//cout<<"inner loop"<<endl;
		//	cout<<"using Lc Lr "<<init_Lc<<" "<<init_Lr<<endl;
			orbit[0] = init_semi;
			orbit[1] = init_Lc;
			orbit[2] = init_Lr;
			forwardSolver(stellar_system_diff_eq, stellar_system_jacobian,
					start_age, orbit, stop_evolution, (void*)this);

	//		cout<<"last age "<<forwardSolver.get_tabulated_indep_var()->back();
			assert(forwardSolver.get_tabulated_indep_var()->back() == curr_age);
			double diff_Lc = curr_spin*curr_Ic -
					forwardSolver.get_tabulated_dep_var(1)->back();
			if (abs(diff_Lc) < 1e-8) break;
			cout<<"wanted actual " << curr_spin*curr_Ic << " "<<
					forwardSolver.get_tabulated_dep_var(1)->back() <<endl;
			init_Lc += diff_Lc;
//			cout<<" diff Lc "<<diff_Lc<<endl;
			init_Lr = init_Lc*init_Ir/init_Ic;
			if (init_Lc <= 0) break;
		}
		valarray<double> new_ages = list_to_valarray(
				*(forwardSolver.get_tabulated_indep_var()));
		valarray<double> new_Lc = list_to_valarray(
				*(forwardSolver.get_tabulated_dep_var(1)));
		valarray<double> new_Lr = list_to_valarray(
				*(forwardSolver.get_tabulated_dep_var(2)));
		valarray<double> new_Lc_deriv = list_to_valarray(
				*(forwardSolver.get_tabulated_dep_var_deriv(1)));
		valarray<double> new_Lr_deriv = list_to_valarray(
				*(forwardSolver.get_tabulated_dep_var_deriv(2)));
		star->set_angular_momentum_evolution(
				new_ages, new_Lc, new_Lc_deriv, convective);
		star->set_angular_momentum_evolution(
				new_ages, new_Lr, new_Lr_deriv, radiative);
		//star->set_angular_momentum_evolution()
	}
}*/

void StellarSystem::output_evolution(std::string &) const
{
	throw Error::Runtime("Method output_evolution of StellarSystem is not implemented yet.");
}

double no_planet_dwconv_dt(double age, const std::valarray<double> &orbit,
		const StellarSystem &system)
{
	double no_planet_diff_eq[2], no_planet_orbit[2];
	double semimajor=orbit[0]*Rsun_AU,
		   worb=system.get_planet().orbital_angular_velocity_semimajor(
				   semimajor),
		   Iconv=system.get_star().moment_of_inertia(age, convective),
		   dIconv_dt=system.get_star().moment_of_inertia_deriv(age,
				   convective);
	no_planet_orbit[0]=worb*Iconv;
	no_planet_orbit[1]=orbit[1];
	system.no_planet_differential_equation(age, no_planet_orbit,
			no_planet_diff_eq);
	return (no_planet_diff_eq[0] - dIconv_dt*worb)/Iconv;
}
