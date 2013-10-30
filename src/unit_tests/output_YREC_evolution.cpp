/**\file
 *
 * \brief Outputs the interpolated stellar quantities derived from the YREC
 * tracks.
 *
 * \ingroup UnitTests_group
 */

#include "output_YREC_evolution.h"

///\brief Outputs the function value, its first and second derivative as
///consecutive real value to os.
void output_all_deriv(std::ostream &os, const EvolvingStellarQuantity *q,
		double age)
{
	const FunctionDerivatives *deriv=q->deriv(age);
	os << std::setw(25) << deriv->order(0)
		<< std::setw(25) << deriv->order(1)
		<< std::setw(25) << deriv->order(2);
	delete deriv;
}

void output_evolution(std::ostream &os,
		const StellarEvolution &evolution, double mass)
{
	std::ios_base::fmtflags orig_os_flags=os.flags();
	std::streamsize orig_os_width=os.width();
	char orig_os_fill=os.fill();
	os.precision(16);
	os.setf(std::ios_base::scientific);
	os.fill(' ');
	os << "#" << std::setw(24) << "age"
		<< std::setw(25) << "R*"
		<< std::setw(25) << "R*_deriv"
		<< std::setw(25) << "R*_2deriv"
		<< std::setw(25) << "L*"
		<< std::setw(25) << "L*_deriv"
		<< std::setw(25) << "L*_2deriv"
		<< std::setw(25) << "Iconv"
		<< std::setw(25) << "Iconv_deriv"
		<< std::setw(25) << "Iconv_2deriv"
		<< std::setw(25) << "Itot"
		<< std::setw(25) << "Itot_deriv"
		<< std::setw(25) << "Itot_2deriv"
		<< std::setw(25) << "Mrad"
		<< std::setw(25) << "Mrad_deriv"
		<< std::setw(25) << "Mrad_2deriv"
		<< std::setw(25) << "Rcore"
		<< std::setw(25) << "Rcore_deriv"
		<< std::setw(25) << "Rcore_2deriv"
		<< std::setw(25) << "Irad"
		<< std::setw(25) << "Irad_deriv"
		<< std::setw(25) << "Irad_2deriv"
		<< std::endl;
	os << "#";
	for(unsigned col=1; col<=19; col++) {
		std::ostringstream colnum;
		colnum << "[" << col << "]";
		os << std::setw((col==1 ? 24 : 25)) << colnum.str();
	}
	os << std::endl;
	const EvolvingStellarQuantity
		*Rstar=evolution.interpolate_radius(mass),
		*Lstar=evolution.interpolate_luminosity(mass),
		*Iconv=evolution.interpolate_moment_of_inertia(mass, convective),
		*Irad=evolution.interpolate_moment_of_inertia(mass, radiative),
		*Itot=evolution.interpolate_moment_of_inertia(mass, total),
		*Mrad=evolution.interpolate_zone_mass(mass, radiative),
		*Rcore=evolution.interpolate_core_boundary(mass);
	double min_age=Rstar->range_low()*(1.0 +
			std::numeric_limits<double>::epsilon()),
		   max_age=Rstar->range_high()*(1.0 -
			std::numeric_limits<double>::epsilon()),
		   age_step=std::pow(max_age/min_age, 1e-3);
	std::cout << "R range: (" << Rstar->range_low() << ", " 
		<< Rstar->range_high() << ")" << std::endl;
	std::cout << "Iconv range: (" << Iconv->range_low() << ", "
		<< Iconv->range_high() << ")" << std::endl;
	std::cout << "Iconv(1e-5)=" << (*Iconv)(1e-5) << std::endl;
	for(double age=min_age; age<=max_age; age*=age_step) {
		os << std::setw(25) << age;
		output_all_deriv(os, Rstar, age);
		output_all_deriv(os, Lstar, age);
		output_all_deriv(os, Iconv, age);
		output_all_deriv(os, Itot, age);
		output_all_deriv(os, Mrad, age);
		output_all_deriv(os, Rcore, age);
		output_all_deriv(os, Irad, age);
		os << std::endl;
	}
	delete Rstar;
	delete Iconv;
	delete Itot;
	delete Mrad;
	delete Rcore;
	os.flags(orig_os_flags);
	os.width(orig_os_width);
	os.fill(orig_os_fill);
}

///Outputs the interpolated stellar quantities derived from the YREC tracks.
int main(int argc, char **argv)
{
	YRECEvolution toSave("../YREC", 0, 2.0, 2.0);
	toSave.save_state("interp_state_data_phs4");
	return 0;
	assert(argc==2);
	try {
//		YRECEvolution evolution("../YREC", 0, 2.0, 2.0);
		YRECEvolution evolution;
		evolution.load_state("../interp_state_data");
		for(double m=0.5; m<1.22; m+=0.001) {
			std::ostringstream fname;
			fname.precision(3);
			fname.setf(std::ios::fixed,std::ios::floatfield);
			fname << argv[1] << "_M" << m << ".evol"; 
			std::cout << "creating file: " << fname.str() << std::endl;
			std::ofstream outf(fname.str().c_str(), std::ios::out);
			output_evolution(outf, evolution, m);
			outf.close();
		}
		const EvolvingStellarQuantity
			*Iconv=evolution.interpolate_moment_of_inertia(1.012, convective);
		std::cerr << "Iconv(m=1.012, t=8.756Gyrs)=" << (*Iconv)(8.756) << std::endl;
	} catch(Error::General &ex) {
		std::cerr << "output_YREC_evolution::main: unexpected exception "
			"thrown: " << ex.what() << ": " << ex.get_message() << std::endl;
	} catch(std::exception &ex) {
		std::cout << "output_YREC_evolution::main: unexpected exception "
			"thrown: " << ex.what() << std::endl;
	}
}
