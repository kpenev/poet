/**\file
 *
 * \brief Defines some of the methods of the classes used to exercise
 * StellarSystem.
 *
 * \ingroup UnitTests_group
 */

#include "test_StellarSystem.h"

const double
	convertI = AstroConst::solar_mass*std::pow(AstroConst::solar_radius, 2),
	convert_torque=AstroConst::day*AstroConst::Gyr/convertI;


test_StellarSystem::test_StellarSystem()
{
	ssdata = new SystemData();
	//ssdata();
	ssdata->create_random_system(&system);
	//return;
	sdata = ssdata->sdata;
	pdata = ssdata->pdata;
//	TEST_ADD(test_StellarSystem::test_basic);
	TEST_ADD(test_StellarSystem::test_unlocked_orbit);
	TEST_ADD(test_StellarSystem::test_locked_orbit);
	TEST_ADD(test_StellarSystem::test_no_planet_orbit);
	TEST_ADD(test_StellarSystem::test_disk_locked_orbit);
	TEST_ADD(test_StellarSystem::test_unlocked_jacobian);
	TEST_ADD(test_StellarSystem::test_locked_jacobian);
	TEST_ADD(test_StellarSystem::test_no_planet_jacobian);
	TEST_ADD(test_StellarSystem::test_disk_locked_jacobian);
//	TEST_ADD(test_StellarSystem::test_jacobian);
}

void test_StellarSystem::setup()
{
}

void test_StellarSystem::tear_down()
{
}

void test_StellarSystem::test_basic()
{
	std::cout << "Age: " << system->current_age() << std::endl;
	std::cout << "Name: " << system->get_name() << std::endl;
}

void DeriveFromOrbit::init_common(double age, double Lrad, const StellarSystem &system)
{
	__Iconv=system.get_star().moment_of_inertia(age, convective)*convertI;
	__Lrad=Lrad*convertI/AstroConst::day;
	__Mstar=system.get_star().get_mass()*AstroConst::solar_mass;
	__Mplanet=system.get_planet().get_mass()*AstroConst::jupiter_mass;
	__Rstar=system.get_star().get_radius(age)*AstroConst::solar_radius;
	__windK=system.get_star().get_wind_strength()*
		convertI*std::pow(AstroConst::day, 2)/AstroConst::Gyr;
	__wind_wsat=system.get_star().get_wind_saturation_frequency()/
		AstroConst::day;
	__coupling_timescale=system.get_star().get_core_env_coupling_timescale()*
		AstroConst::Gyr;
	__dIconv_dt=system.get_star().moment_of_inertia_deriv(age, convective)*
		convertI/AstroConst::Gyr;
	if(age>system.get_star().core_formation_age()) {
		__Irad=system.get_star().moment_of_inertia(age, radiative)*convertI;
		__Rrad=system.get_star().get_rad_radius(age)*
			AstroConst::solar_radius;
		__Mrad=system.get_star().get_rad_mass(age)*AstroConst::solar_mass;
		__Mrad_deriv=system.get_star().get_rad_mass_deriv(age)*
			AstroConst::solar_mass/AstroConst::Gyr;
	} else __Irad=__Rrad=__Mrad=__Mrad_deriv=0;

}

void DeriveFromOrbit::operator()(double age, double a, double Lconv,
		double Lrad, const StellarSystem &system)
{
	init_common(age, Lrad, system);
	__a=a*AstroConst::solar_radius;
	__Lconv=Lconv*convertI/AstroConst::day;
	__wconv=__Lconv/__Iconv;
	__worb=sqrt(AstroConst::G*(__Mstar+__Mplanet)/std::pow(__a, 3));
	__differential_rotation=(__Iconv*__Lrad - __Irad*__Lconv)/
		(__Iconv+__Irad);
	__wind_torque=-__windK*__wconv*std::pow(
			std::min(__wconv, __wind_wsat), 2)*
		std::pow(__Rstar/AstroConst::solar_radius, 0.5)*
		std::pow(__Mstar/AstroConst::solar_mass, -0.5);
	__coupling_torque=__differential_rotation/__coupling_timescale;
	if(age>=system.get_star().core_formation_age()) 
		__core_growth_torque=2.0/3*__Rrad*__Rrad*__wconv*__Mrad_deriv;
	else __core_growth_torque=0;
}

void DeriveFromLockedOrbit::operator()(double age, double a, double Lrad,
		const StellarSystem &system)
{
	init_common(age, Lrad, system);
	__a=a*AstroConst::solar_radius;
	__worb=sqrt(AstroConst::G*(__Mstar+__Mplanet)/std::pow(__a, 3));
	__wconv=__worb;
	__Lconv=__wconv*__Iconv;
	__differential_rotation=(__Iconv*__Lrad - __Irad*__Lconv)/
		(__Iconv+__Irad);
	__wind_torque=-__windK*__wconv*std::pow(
			std::min(__wconv, __wind_wsat), 2)*
		std::pow(__Rstar/AstroConst::solar_radius, 0.5)*
		std::pow(__Mstar/AstroConst::solar_mass, -0.5);
	__coupling_torque=__differential_rotation/__coupling_timescale;
	if(age>system.get_star().core_formation_age()) 
		__core_growth_torque=2.0/3*__Rrad*__Rrad*__wconv*__Mrad_deriv;
	else __core_growth_torque=0;
}

void DeriveFromNoPlanet::operator()(double age, double Lconv, double Lrad,
		const StellarSystem &system)
{
	init_common(age, Lrad, system);
	__a=NaN;
	__Lconv=Lconv*convertI/AstroConst::day;
	__wconv=__Lconv/__Iconv;
	__worb=NaN;
	__differential_rotation=(__Iconv*__Lrad - __Irad*__Lconv)/
		(__Iconv+__Irad);
	__wind_torque=-__windK*__wconv*std::pow(
			std::min(__wconv, __wind_wsat), 2)*
		std::pow(__Rstar/AstroConst::solar_radius, 0.5)*
		std::pow(__Mstar/AstroConst::solar_mass, -0.5);
	__coupling_torque=__differential_rotation/__coupling_timescale;
	if(age>=system.get_star().core_formation_age()) 
		__core_growth_torque=2.0/3*__Rrad*__Rrad*__wconv*__Mrad_deriv;
	else __core_growth_torque=0;
}

void DeriveFromDiskLocked::operator()(double age, double Lrad,
		const StellarSystem &system)
{
	init_common(age, Lrad, system);
	__a=NaN;
	__wconv=system.get_star().get_disk_lock_frequency()/AstroConst::day;
	__Lconv=__Iconv*__wconv;
	__worb=NaN;
	__differential_rotation=(__Iconv*__Lrad - __Irad*__Lconv)/
		(__Iconv+__Irad);
	__wind_torque=-__windK*__wconv*std::pow(
			std::min(__wconv, __wind_wsat), 2)*
		std::pow(__Rstar/AstroConst::solar_radius, 0.5)*
		std::pow(__Mstar/AstroConst::solar_mass, -0.5);
	__coupling_torque=__differential_rotation/__coupling_timescale;
	if(age>=system.get_star().core_formation_age()) 
		__core_growth_torque=2.0/3*__Rrad*__Rrad*__wconv*__Mrad_deriv;
	else __core_growth_torque=0;
}

/*
void test_StellarSystem::predict_result(int ageIndex, long double* results)
{
	//results: a^6.5, Lconv, Lrad
	using namespace AstroConst;
	long double age = sdata->evolution_ages[ageIndex]*Gyr;
	long double Lconv = sdata->Lconv[ageIndex]*convertI/day;
	long double Lrad = sdata->Lrad[ageIndex]*convertI/day;
	long double Iconv = sdata->Iconv[ageIndex]*convertI;
	long double Irad = sdata->Irad[ageIndex]*convertI;

	long double delta_L = (Iconv*Lrad-Irad*Lconv)/(Iconv+Irad);
	long double timescale = sdata->coupling_timescale*Gyr;

	long double wconv = Lconv/Iconv;
	long double wsat = sdata->wind_sat_freq/day;
	long double semi = pdata->get_semi(age/Gyr)*AU;
	long double Rrad = sdata->Rrad[ageIndex]*solar_radius;

	long double pmass = pdata->mass*jupiter_mass;
	long double smass = sdata->mass*solar_mass;

	long double sradius = sdata->all_radii[ageIndex]*solar_radius;
	long double worb = sqrt(G*smass/std::pow(semi, 3));

	long double Mrad_deriv = sdata->Mrad_deriv[ageIndex]*solar_mass/Gyr;
	long double sign;
	if (wconv > worb) sign = 1;
	else sign = -1;
	double Q;
	double tidal_freq = (wconv - worb)*day;
	if (abs(tidal_freq) < sdata->Q_trans_width)
		Q = sdata->tidal_Q/sin(M_PI/2.0/sdata->Q_trans_width*tidal_freq);
	else Q = sdata->tidal_Q;
	results[0] = sign*4.5*sqrt(G/smass)*std::pow(sradius, 5)*pmass/
		std::pow(semi, (long double)5.5)/Q;
	long double dLc_tide = -0.5*pmass*smass*sqrt(G/semi/(smass+pmass))*
		results[0];
	long double K = sdata->wind_strength*convertI*day*day/Gyr;
	long double dLc_wind = -K*wconv*std::pow(min(wconv,wsat), 2)*
		std::pow(sradius/solar_radius, (long double)0.5)*
		std::pow(smass/solar_mass, (long double)-0.5);

	long double mutual_torque = delta_L/timescale - 2.0/3*Rrad*Rrad*wconv*Mrad_deriv;
	//		long double core_gain = 2.0/3*Rrad*Rrad*Mrad_deriv;

	long double dLc_tot = mutual_torque + dLc_tide + dLc_wind;
	long double dLr_tot = -mutual_torque;
	results[1] = dLc_tot;
	results[2] = dLr_tot;

	//now convert to right units
	results[0] *= 6.5*std::pow(semi, (long double)5.5);
	results[0] /= std::pow(solar_radius,6.5)/Gyr;
	results[1] /= convertI/day/Gyr;
	results[2] /= convertI/day/Gyr;
}*/

void test_StellarSystem::predict_unlocked_orbit_diff_eq(
		const DeriveFromOrbit &quantities_SI, double* test_derivs,
		int dadt_sign)
{
	double tidal_freq = (quantities_SI.wconv() - quantities_SI.worb())*
		AstroConst::day, Q;
	int sign;
	if(dadt_sign==0) {
		sign=(tidal_freq>0 ? 1 : -1);
		if (abs(tidal_freq) < sdata->Q_trans_width)
			Q = sdata->tidal_Q/sin(M_PI/2.0/sdata->Q_trans_width*tidal_freq);
		else Q = sdata->tidal_Q;
	} else {
		sign=dadt_sign;
		Q=sdata->tidal_Q;
	}

	//Equations in SI
	test_derivs[0]=4.5*sign*sqrt(AstroConst::G/quantities_SI.Mstar())*
		std::pow(quantities_SI.Rstar(), 5)*quantities_SI.Mplanet()/
		std::pow(quantities_SI.a(), 5.5)/Q;
	double tidal_torque = -0.5*quantities_SI.Mplanet()*
		quantities_SI.Mstar()*sqrt(AstroConst::G/quantities_SI.a()/
				(quantities_SI.Mstar()+quantities_SI.Mplanet()))*
		test_derivs[0];

	test_derivs[2]=quantities_SI.core_growth_torque()-
		quantities_SI.coupling_torque();
	test_derivs[1]=quantities_SI.wind_torque() + tidal_torque -
		test_derivs[2];

	//Convert to orbital evolution units
	test_derivs[0] *= 6.5*std::pow(quantities_SI.a(), 5.5)/
		std::pow(AstroConst::solar_radius,6.5)*AstroConst::Gyr;
	test_derivs[1] *= convert_torque;
	test_derivs[2] *= convert_torque;
}

void test_StellarSystem::predict_locked_orbit_diff_eq(
		const DeriveFromOrbit &quantities_SI, double* test_derivs)
{
	//Equations in SI
	test_derivs[0]=2.0*(quantities_SI.coupling_torque() -
			quantities_SI.core_growth_torque() +
			quantities_SI.wind_torque() -
			quantities_SI.wconv()*quantities_SI.dIconv_dt())/
		(quantities_SI.Mstar()*quantities_SI.Mplanet()*
		 std::sqrt(AstroConst::G/(quantities_SI.a()*
				 (quantities_SI.Mstar()+quantities_SI.Mplanet()))) -
		 3.0*quantities_SI.Iconv()*quantities_SI.wconv()/quantities_SI.a());
	test_derivs[1]=quantities_SI.core_growth_torque()-
		quantities_SI.coupling_torque();

	//Convert to orbital evolution units
	test_derivs[0] *= AstroConst::Gyr/AstroConst::solar_radius;
	test_derivs[1] *= convert_torque;
}

void test_StellarSystem::predict_no_planet_diff_eq(
		const DeriveFromOrbit &quantities_SI, double *test_derivs)
{
	//Equations in SI
	test_derivs[1]=quantities_SI.core_growth_torque()-
		quantities_SI.coupling_torque();
	test_derivs[0]=quantities_SI.wind_torque()-test_derivs[1];

	//Convert to orbital evolution units
	test_derivs[0] *= convert_torque;
	test_derivs[1] *= convert_torque;
}

void test_StellarSystem::predict_disk_locked_diff_eq(
		const DeriveFromOrbit &quantities_SI, double *test_derivs)
{
	//Equation in SI
	test_derivs[0] = quantities_SI.core_growth_torque()-
		quantities_SI.coupling_torque();

	//Convert to orbital evolution units
	test_derivs[0] *= convert_torque;
}

void test_StellarSystem::test_interpolate()
{
	TEST_ASSERT(sdata->evolution_ages.size() == sdata->all_radii.size());
	for (size_t i=0; i < sdata->evolution_ages.size(); i++) {
		double age = sdata->evolution_ages[i];
		double interp_radius = ssdata->star->get_radius(age);
		double real_radius = sdata->all_radii[i];
		std::cout << "age " << age << ": " << interp_radius << " " 
			<< real_radius << std::endl;
	}
}

void test_StellarSystem::test_unlocked_orbit()
{
	for (size_t i=0; i < sdata->evolution_ages.size(); i++) {
		double age = sdata->evolution_ages[i];
		if (age >= system->get_planet().get_lifetime()) break;
		double Lconv = sdata->Lconv[i],
			   Lrad = (age>=system->get_star().core_formation_age() ?
					   sdata->Lrad[i] : 0),
			   a=pdata->get_semi(age)*AstroConst::AU/AstroConst::solar_radius,
			   a6p5 = std::pow(a, 6.5),
			   params[3], derivs[3], test_derivs[3];
		params[0]=a6p5; params[1]=Lconv; params[2]=Lrad;
		derivs[0]=0; derivs[1]=0; derivs[2]=0;

		DeriveFromOrbit quantities_SI(age, a, Lconv, Lrad, *system);
		for(int sign=-1; sign<=1; sign++) {
			system->orbit_differential_equation(age, params, derivs, sign);
			predict_unlocked_orbit_diff_eq(quantities_SI, test_derivs, sign);

			std::ostringstream msg;

			msg << "StellarSystem unlocked orbit derivatives for age=" << age
				<< ", sign=" << sign << ": " << derivs[0] << ", "
				<< derivs[1] << ", " << derivs[2] << ", manual derivatives: "
				<< test_derivs[0] << ", " << test_derivs[1] << ", "
				<< test_derivs[2];
			TEST_ASSERT_MSG(check_diff(derivs[0], test_derivs[0], 1e-12, 0)&&
					check_diff(derivs[1], test_derivs[1], 1e-12, 0) &&
					check_diff(derivs[2], test_derivs[2], 1e-12, 0),
					msg.str().c_str())
		}
	}
}

void test_StellarSystem::test_locked_orbit()
{
	try {
		for (size_t i=0; i < sdata->evolution_ages.size(); i++) {
			double age = sdata->evolution_ages[i];
			if (age >= system->get_planet().get_lifetime()) break;
			double Lrad = (age>=system->get_star().core_formation_age() ?
					sdata->Lrad[i] : 0),
				   a=pdata->get_semi(age)*AstroConst::AU/
					   AstroConst::solar_radius,
				   params[2], derivs[2], test_derivs[2];
			params[0]=a; params[1]=Lrad;
			derivs[0]=0; derivs[1]=0;
			system->locked_orbit_differential_equation(age, params,
					derivs);
			DeriveFromLockedOrbit quantities_SI(age, a, Lrad, *system);
			predict_locked_orbit_diff_eq(quantities_SI, test_derivs);

			std::ostringstream msg;
			msg << "StellarSystem derivatives: " << derivs[0] << ", "
				<< derivs[1]  << ", manual derivatives: "
				<< test_derivs[0] << ", " <<test_derivs[1];
			TEST_ASSERT_MSG(
					check_diff(derivs[0], test_derivs[0], 1e-12, 0) &&
					check_diff(derivs[1], test_derivs[1], 1e-12, 0),
					msg.str().c_str())
		}
	} catch (const std::exception &ex) {
		std::ostringstream msg;
		msg << "Unexpected exception thrown in test_locked_orbit: "
			<< ex.what();
		TEST_ASSERT_MSG(false, msg.str().c_str());
	}
}

void test_StellarSystem::test_no_planet_orbit()
{
	try {
		for (size_t i=0; i < sdata->evolution_ages.size(); i++) {
			double age = sdata->evolution_ages[i];
			if (age >= system->get_planet().get_lifetime()) break;
			double Lconv = sdata->Lconv[i],
				   Lrad = (age>=system->get_star().core_formation_age() ?
						   sdata->Lrad[i] : 0),
				   params[2], derivs[2], test_derivs[2];
			params[0]=Lconv; params[1]=Lrad;
			derivs[0]=0; derivs[1]=0;
			system->no_planet_differential_equation(age, params,
					derivs);
			DeriveFromNoPlanet quantities_SI(age, Lconv, Lrad, *system);
			predict_no_planet_diff_eq(quantities_SI, test_derivs);

			std::ostringstream msg;
			msg << "StellarSystem derivatives: " << derivs[0] << ", "
				<< derivs[1]  << ", manual derivatives: "
				<< test_derivs[0] << ", " <<test_derivs[1];
			TEST_ASSERT_MSG(
					check_diff(derivs[0], test_derivs[0], 1e-12, 0) &&
					check_diff(derivs[1], test_derivs[1], 1e-12, 0),
					msg.str().c_str())
		}
	} catch (const std::exception &ex) {
		std::ostringstream msg;
		msg << "Unexpected exception thrown in test_locked_orbit: "
			<< ex.what();
		TEST_ASSERT_MSG(false, msg.str().c_str());
	}
}

void test_StellarSystem::test_disk_locked_orbit()
{
	try {
		for (size_t i=0; i < sdata->evolution_ages.size(); i++) {
			double age = sdata->evolution_ages[i];
			if (age >= system->get_planet().get_lifetime()) break;
			double Lrad = (age>=system->get_star().core_formation_age() ?
						   sdata->Lrad[i] : 0),
				   params[1], derivs[1], test_derivs[1];
			params[0]=Lrad;
			derivs[0]=0;
			system->locked_conv_differential_equation(age, params,
					derivs);
			DeriveFromDiskLocked quantities_SI(age, Lrad, *system);
			predict_disk_locked_diff_eq(quantities_SI, test_derivs);

			std::ostringstream msg;
			msg << "StellarSystem derivatives: " << derivs[0]
				<< ", manual derivatives: " << test_derivs[0];
			TEST_ASSERT_MSG(check_diff(derivs[0], test_derivs[0], 1e-12, 0),
					msg.str().c_str())
		}
	} catch (const std::exception &ex) {
		std::ostringstream msg;
		msg << "Unexpected exception thrown in test_locked_orbit: "
			<< ex.what();
		TEST_ASSERT_MSG(false, msg.str().c_str());
	}
}

void test_StellarSystem::test_unlocked_jacobian()
{
	try {
		for (size_t i=1; i < sdata->evolution_ages.size(); i++) {
			double age = sdata->evolution_ages[i];
			if (age >= system->get_planet().get_lifetime()) break;
			double Lconv = sdata->Lconv[i],
				   Lrad = (age>=system->get_star().core_formation_age() ?
						   sdata->Lrad[i] : 0),
				   a=pdata->get_semi(age)*
					   AstroConst::AU/AstroConst::solar_radius,
				   a6p5 = std::pow(a, 6.5),
				   age_step=1e-6;
			std::valarray<double> params(3), derivs_plus(3), derivs_minus(3),
				param_steps(1e-4, 3), jacobian(9), estimated_jacobian(9),
				time_derivs(3), estimated_time_derivs(3);
			params[0]=a6p5; params[1]=Lconv; params[2]=Lrad;

			for(int sign=-1; sign<=1; sign++) {
				for(size_t param_ind=0; param_ind<params.size();
						param_ind++) {
					std::valarray<double> params_minus=params,
						params_plus=params;
					params_minus[param_ind]-=param_steps[param_ind];
					params_plus[param_ind]+=param_steps[param_ind];
					system->orbit_differential_equation(age,
							&(params_plus[0]), &(derivs_plus[0]), sign);
					system->orbit_differential_equation(age,
							&(params_minus[0]), &(derivs_minus[0]), sign);
					estimated_jacobian[std::slice(param_ind,
							params.size(), params.size())]=
						0.5*(derivs_plus-derivs_minus)/
						param_steps[param_ind];
				}
				system->orbit_differential_equation(age+age_step,
						&(params[0]), &(derivs_plus[0]), sign);
				system->orbit_differential_equation(age-age_step,
						&(params[0]), &(derivs_minus[0]), sign);
				estimated_time_derivs=
					0.5*(derivs_plus-derivs_minus)/age_step;

				system->orbit_jacobian(age, &(params[0]), &(jacobian[0]),
						&(time_derivs[0]), sign);

				std::ostringstream msg;

				msg << "StellarSystem unlocked orbit jacobian for age="
					<< age << ", sign=" << sign << ": " << jacobian
					<< ", estimated jacobian: " << estimated_jacobian;
				TEST_ASSERT_MSG(
						check_diff(jacobian, estimated_jacobian,
							std::valarray<double>(1e-5, jacobian.size()),
							std::valarray<double>(1e-8, jacobian.size())),
						msg.str().c_str());

				msg.str("");
				msg << "StellarSystem unlocked orbit time derivatives for "
				"age=" << age << ", sign=" << sign << ": " << time_derivs
					<< ", estimated time derivatives: "
					<< estimated_time_derivs;
				TEST_ASSERT_MSG(
						check_diff(time_derivs, estimated_time_derivs,
							std::valarray<double>(1e-6, time_derivs.size()),
							std::valarray<double>(1e-8, time_derivs.size())),
						msg.str().c_str())
			}
		}
	} catch (const std::exception &ex) {
		std::ostringstream msg;
		msg << "Unexpected exception thrown in test_unlocked_jacobian: "
			<< ex.what();
		TEST_ASSERT_MSG(false, msg.str().c_str());
	}
}

void test_StellarSystem::test_locked_jacobian()
{
	try {
		for (size_t i=1; i < sdata->evolution_ages.size(); i++) {
			double age = sdata->evolution_ages[i];
			if (age >= system->get_planet().get_lifetime()) break;
			double Lrad = (age>=system->get_star().core_formation_age() ?
						   sdata->Lrad[i] : 0),
				   a=pdata->get_semi(age)*
					   AstroConst::AU/AstroConst::solar_radius,
				   age_step=1e-6;
			std::valarray<double> params(2), derivs_plus(2), derivs_minus(2),
				param_steps(1e-4, 2), jacobian(4), estimated_jacobian(4),
				time_derivs(2), estimated_time_derivs(2);
			params[0]=a; params[1]=Lrad;

			for(size_t param_ind=0; param_ind<params.size();
					param_ind++) {
				std::valarray<double> params_minus=params,
					params_plus=params;
				params_minus[param_ind]-=param_steps[param_ind];
				params_plus[param_ind]+=param_steps[param_ind];
				system->locked_orbit_differential_equation(age,
						&(params_plus[0]), &(derivs_plus[0]));
				system->locked_orbit_differential_equation(age,
						&(params_minus[0]), &(derivs_minus[0]));
				estimated_jacobian[std::slice(param_ind,
						params.size(), params.size())]=
					0.5*(derivs_plus-derivs_minus)/
					param_steps[param_ind];
			}
			system->locked_orbit_differential_equation(age+age_step,
					&(params[0]), &(derivs_plus[0]));
			system->locked_orbit_differential_equation(age-age_step,
					&(params[0]), &(derivs_minus[0]));
			estimated_time_derivs=
				0.5*(derivs_plus-derivs_minus)/age_step;

			system->locked_orbit_jacobian(age, &(params[0]), &(jacobian[0]),
					&(time_derivs[0]));

			std::ostringstream msg;

			msg << "StellarSystem locked orbit jacobian for age="
				<< age << ": " << jacobian
				<< ", estimated jacobian: " << estimated_jacobian;
			TEST_ASSERT_MSG(
					check_diff(jacobian, estimated_jacobian,
						std::valarray<double>(1e-5, jacobian.size()),
						std::valarray<double>(1e-8, jacobian.size())),
					msg.str().c_str());

			msg.str("");
			msg << "StellarSystem locked orbit time derivatives for "
				"age=" << age << ": " << time_derivs
				<< ", estimated time derivatives: "
				<< estimated_time_derivs;
			TEST_ASSERT_MSG(
					check_diff(time_derivs, estimated_time_derivs,
						std::valarray<double>(1e-6, time_derivs.size()),
						std::valarray<double>(1e-8, time_derivs.size())),
					msg.str().c_str())
		}
	} catch (const std::exception &ex) {
		std::ostringstream msg;
		msg << "Unexpected exception thrown in test_unlocked_jacobian: "
			<< ex.what();
		TEST_ASSERT_MSG(false, msg.str().c_str());
	}
}

void test_StellarSystem::test_no_planet_jacobian()
{
	try {
		for (size_t i=1; i < sdata->evolution_ages.size(); i++) {
			double age = sdata->evolution_ages[i];
			if (age >= system->get_planet().get_lifetime()) break;
			double Lconv=sdata->Lconv[i],
				   Lrad = (age>=system->get_star().core_formation_age() ?
						   sdata->Lrad[i] : 0),
				   age_step=1e-6;
			std::valarray<double> params(2), derivs_plus(2), derivs_minus(2),
				param_steps(1e-4, 2), jacobian(4), estimated_jacobian(4),
				time_derivs(2), estimated_time_derivs(2);
			params[0]=Lconv; params[1]=Lrad;

			for(size_t param_ind=0; param_ind<params.size();
					param_ind++) {
				std::valarray<double> params_minus=params,
					params_plus=params;
				params_minus[param_ind]-=param_steps[param_ind];
				params_plus[param_ind]+=param_steps[param_ind];
				system->no_planet_differential_equation(age,
						&(params_plus[0]), &(derivs_plus[0]));
				system->no_planet_differential_equation(age,
						&(params_minus[0]), &(derivs_minus[0]));
				estimated_jacobian[std::slice(param_ind,
						params.size(), params.size())]=
					0.5*(derivs_plus-derivs_minus)/
					param_steps[param_ind];
			}
			system->no_planet_differential_equation(age+age_step,
					&(params[0]), &(derivs_plus[0]));
			system->no_planet_differential_equation(age-age_step,
					&(params[0]), &(derivs_minus[0]));
			estimated_time_derivs=
				0.5*(derivs_plus-derivs_minus)/age_step;

			system->no_planet_jacobian(age, &(params[0]), &(jacobian[0]),
					&(time_derivs[0]));

			std::ostringstream msg;

			msg << "StellarSystem no planet jacobian for age="
				<< age << ": " << jacobian
				<< ", estimated jacobian: " << estimated_jacobian;
			TEST_ASSERT_MSG(
					check_diff(jacobian, estimated_jacobian,
						std::valarray<double>(1e-5, jacobian.size()),
						std::valarray<double>(1e-8, jacobian.size())),
					msg.str().c_str());

			msg.str("");
			msg << "StellarSystem no planet time derivatives for "
				"age=" << age << ": " << time_derivs
				<< ", estimated time derivatives: "
				<< estimated_time_derivs;
			TEST_ASSERT_MSG(
					check_diff(time_derivs, estimated_time_derivs,
						std::valarray<double>(1e-6, time_derivs.size()),
						std::valarray<double>(1e-8, time_derivs.size())),
					msg.str().c_str())
		}
	} catch (const std::exception &ex) {
		std::ostringstream msg;
		msg << "Unexpected exception thrown in test_unlocked_jacobian: "
			<< ex.what();
		TEST_ASSERT_MSG(false, msg.str().c_str());
	}
}

void test_StellarSystem::test_disk_locked_jacobian()
{
	try {
		for (size_t i=1; i < sdata->evolution_ages.size(); i++) {
			double age = sdata->evolution_ages[i];
			if (age >= system->get_planet().get_lifetime()) break;
			double Lrad = (age>=system->get_star().core_formation_age() ?
						   sdata->Lrad[i] : 0),
				   age_step=1e-6;
			std::valarray<double> params(1), derivs_plus(1), derivs_minus(1),
				param_steps(1e-4, 1), jacobian(1), estimated_jacobian(1),
				time_derivs(1), estimated_time_derivs(1);
			params[0]=Lrad;

			for(size_t param_ind=0; param_ind<params.size();
					param_ind++) {
				std::valarray<double> params_minus=params,
					params_plus=params;
				params_minus[param_ind]-=param_steps[param_ind];
				params_plus[param_ind]+=param_steps[param_ind];
				system->locked_conv_differential_equation(age,
						&(params_plus[0]), &(derivs_plus[0]));
				system->locked_conv_differential_equation(age,
						&(params_minus[0]), &(derivs_minus[0]));
				estimated_jacobian[std::slice(param_ind,
						params.size(), params.size())]=
					0.5*(derivs_plus-derivs_minus)/
					param_steps[param_ind];
			}
			system->locked_conv_differential_equation(age+age_step,
					&(params[0]), &(derivs_plus[0]));
			system->locked_conv_differential_equation(age-age_step,
					&(params[0]), &(derivs_minus[0]));
			estimated_time_derivs=
				0.5*(derivs_plus-derivs_minus)/age_step;

			system->locked_conv_jacobian(age, &(params[0]), &(jacobian[0]),
					&(time_derivs[0]));

			std::ostringstream msg;

			msg << "StellarSystem no planet jacobian for age="
				<< age << ": " << jacobian
				<< ", estimated jacobian: " << estimated_jacobian;
			TEST_ASSERT_MSG(
					check_diff(jacobian, estimated_jacobian,
						std::valarray<double>(1e-5, jacobian.size()),
						std::valarray<double>(1e-8, jacobian.size())),
					msg.str().c_str());

			msg.str("");
			msg << "StellarSystem no planet time derivatives for "
				"age=" << age << ": " << time_derivs
				<< ", estimated time derivatives: "
				<< estimated_time_derivs;
			TEST_ASSERT_MSG(
					check_diff(time_derivs, estimated_time_derivs,
						std::valarray<double>(1e-6, time_derivs.size()),
						std::valarray<double>(1e-8, time_derivs.size())),
					msg.str().c_str())
		}
	} catch (const std::exception &ex) {
		std::ostringstream msg;
		msg << "Unexpected exception thrown in test_unlocked_jacobian: "
			<< ex.what();
		TEST_ASSERT_MSG(false, msg.str().c_str());
	}
}

/*
void test_StellarSystem::test_equations()
{
	using namespace std;
	using namespace AstroConst;
	for (size_t i=0; i < sdata->evolution_ages.size(); i++) {
		double age = sdata->evolution_ages[i];
		if (age >= ssdata->planet->get_lifetime()) break;
		double Lc = sdata->Lconv[i];
		double Lr = sdata->Lrad[i];
		double a6p5 = std::pow(pdata->get_semi(age)*AU/solar_radius,6.5);
		double params[] = {a6p5, Lc, Lr};
		double derivs[3];
		long double testDerivs[3];
		derivs[0] = derivs[1] = derivs[2] = 0;
		system->orbit_differential_equation(age, params, derivs);

		cout << "Real deriv " << derivs[0] << " " <<derivs[1] <<" "
			<<derivs[2]<<endl;
		predict_result(i, testDerivs);
		cout << "Test deriv " << testDerivs[0] << " " <<testDerivs[1] <<
			" "<<testDerivs[2]<<endl;
		cout << endl;
	}
}*/

void test_StellarSystem::check_param_derivs(double age, double* params,
		int index, double* jacobian)
{
	//index is index of varying parameter
	using namespace std;
	double step = params[index]*1e-3;
	double upper[3], lower[3], upperResult[3], lowerResult[3];
	for (int i=0; i < 3; i++) {
		upper[i] = params[i];
		lower[i] = params[i];
	}
	upper[index] += step/2.0;
	lower[index] -= step/2.0;
	system->orbit_differential_equation(age, upper, upperResult);
	system->orbit_differential_equation(age, lower, lowerResult);
	for (int i=0; i < 3; i++) {
		jacobian[3*i + index] = (upperResult[i] - lowerResult[i])/step;
	}
}

void test_StellarSystem::eval_jacobian(double age, double* params,
		double* jacobian, long double* ageDeriv)
{
	using namespace std;
	long double ageStep = age*1e-4;

	double upperResult[3], lowerResult[3];

	system->orbit_differential_equation(age+ageStep/2, params, upperResult);
	system->orbit_differential_equation(age-ageStep/2, params, lowerResult);
	for (int i=0; i < 3; i++)
		ageDeriv[i] = (upperResult[i] - lowerResult[i])/ageStep;
	for (int i=0; i < 3; i++) {
		check_param_derivs(age, params, i, jacobian);
	}
}

void test_StellarSystem::test_jacobian()
{
	using namespace std;
	using namespace AstroConst;
	cout<<sdata->evolution_ages.size()<<endl;
	for (size_t i=1; i < sdata->evolution_ages.size(); i++) {
		double age = sdata->evolution_ages[i];
		if (age >= ssdata->planet->get_lifetime()) break;
		double a6p5 = std::pow(pdata->get_semi(age)*AU/solar_radius,6.5);
		double Lc = sdata->Lconv[i];
		double Lr = sdata->Lrad[i];
		double params[] = {a6p5, Lc, Lr};
		double ageDeriv[3];
		long double testAgeDeriv[3];
		double paramDerivs[9];
		double testParamDerivs[9];

		system->orbit_jacobian(age, params, paramDerivs, ageDeriv);
		eval_jacobian(age, params, testParamDerivs, testAgeDeriv);
		cout << "Real Jacobian age derivs " <<ageDeriv[0]<<" "
			<<ageDeriv[1]<<" "<<ageDeriv[2]<<endl;
		cout<<"Test Jacobian age derivs "<<testAgeDeriv[0]<<" "<<
			testAgeDeriv[1]<<" "<<testAgeDeriv[2]<<endl<<endl;
		cout << "Real Jacobian param derivs ";
		for (int j=0; j < 9; j++) {
			cout << paramDerivs[j]<<" ";
		}
		cout<<endl<<"Test Jacobian param derivs ";
		for (int j=0; j < 9; j++) {
			cout << testParamDerivs[j]<<" ";
		}
		cout<<endl<<endl;
	}
}

#ifdef STANDALONE
int main()
{
	Test::TextOutput output(Test::TextOutput::Verbose);
	test_StellarSystem tests;
	return (tests.run(output) ? EXIT_SUCCESS : EXIT_FAILURE);

}
#endif
