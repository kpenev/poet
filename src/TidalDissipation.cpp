/**\file
 *
 * \brief Declarations of some methods of the TidalDissipation class.
 *
 * \ingroup StellarSystem_group
 */

#include "TidalDissipation.h"

const double TidalDissipation::__Umm_coef[][3]={
	{std::sqrt(3.0*M_PI/10.0)/4.0,
	 -std::sqrt(6.0*M_PI/5.0)/4.0,
	 std::sqrt(3.0*M_PI/10.0)/4.0},

	{-std::sqrt(3.0*M_PI/10.0)/2.0,
	 std::sqrt(6.0*M_PI/5.0)/2.0,
	 std::sqrt(3.0*M_PI/10.0)/2.0},

	{3.0*std::sqrt(M_PI/5.0)/4.0,
	 -std::sqrt(M_PI/5.0)/2.0,
	 3.0*std::sqrt(M_PI/5.0)/4.0},

	{-std::sqrt(3.0*M_PI/10.0)/2.0,
	 std::sqrt(6.0*M_PI/5.0)/2.0,
	 std::sqrt(3.0*M_PI/10.0)/2.0},
									   
	{std::sqrt(3.0*M_PI/10.0)/4.0,
	 -std::sqrt(6.0*M_PI/5.0)/4.0,
	 std::sqrt(3.0*M_PI/10.0)/4.0}};

const double TidalDissipation::__torque_x_plus_coef[]={1.0, 
													   std::sqrt(1.5),
													   std::sqrt(1.5),
													   1.0,
													   0.0};

const double TidalDissipation::__torque_x_minus_coef[]={0.0, 
													    1.0,
													    std::sqrt(1.5),
													    std::sqrt(1.5),
													    1.0};

std::ostream &operator<<(std::ostream &os, 
		const Dissipation::Quantity &quantity)
{
	switch(quantity) {
		case Dissipation::POWER : os << "POWER"; break;
		case Dissipation::TORQUEX : os << "TORQUEX"; break;
		case Dissipation::TORQUEZ : os << "TORQUEZ"; break;
		case Dissipation::SEMIMAJOR_DECAY : os << "SEMIMAJOR_DECAY"; break;
		case Dissipation::ORBIT_SPINUP : os << "ORBIT_SPINUP"; break;
		case Dissipation::INCLINATION_DECAY : os << "INCLINATION_DECAY";
											  break;
		default :
#ifdef DEBUG
				  assert(false)
#endif
					  ;
	};
	return os;
}

///More civilized output for Dissipation::Derivative variables.
std::ostream &operator<<(std::ostream &os,
		const Dissipation::Derivative &deriv)
{
	switch(deriv) {
		case Dissipation::NO_DERIV : std::cout << "NO_DERIV"; break;
		case Dissipation::AGE : std::cout << "AGE"; break;
		case Dissipation::RADIUS : std::cout << "RADIUS"; break;
		case Dissipation::MOMENT_OF_INERTIA :
								   std::cout << "MOMENT_OF_INERTIA";
								   break;
		case Dissipation::SPIN_ANGMOM : std::cout << "SPIN_ANGMOM"; break;
		case Dissipation::SEMIMAJOR : std::cout << "SEMIMAJOR"; break;
		case Dissipation::INCLINATION : std::cout << "INCLINATION"; break;
		default :
#ifdef DEBUG
				  assert(false)
#endif
					  ;
	};
	return os;
}

void TidalDissipation::fill_Umm(double inclination, bool deriv)
{
	double c=std::cos(inclination), s=std::sin(inclination),
		   s2=std::pow(s, 2), sc=s*c, cp1=c+1.0, cm1=c-1.0;

	__Umm[0][0]=__Umm_coef[0][0]*(deriv ? -2.0*s*cp1 : std::pow(cp1, 2));
	__Umm[1][0]=__Umm_coef[1][0]*(deriv ? cp1+2.0*s2 : s*cp1);
	__Umm[2][0]=__Umm_coef[2][0]*(deriv ? 2.0*sc : s2);
	__Umm[3][0]=-__Umm_coef[3][0]*(deriv ? c*cm1-s2 : s*cm1);
	__Umm[4][0]=__Umm_coef[4][0]*(deriv ? -2.0*s*cm1 : std::pow(cm1, 2));

	__Umm[0][1]=__Umm_coef[0][1]*(deriv ? 2.0*sc : s2);
	__Umm[1][1]=__Umm_coef[1][1]*(deriv ? 1.0-2.0*s2 : sc);
	__Umm[2][1]=__Umm_coef[2][1]*(deriv ? -6.0*sc : 2.0-3.0*s2);
	__Umm[3][1]=__Umm_coef[3][1]*(deriv ? 1.0-2.0*s2 : sc);
	__Umm[4][1]=__Umm_coef[4][1]*(deriv ? 2.0*sc : s2);

	__Umm[0][2]=__Umm_coef[0][2]*(deriv ? -2.0*cm1*s : std::pow(cm1, 2));
	__Umm[1][2]=-__Umm_coef[1][2]*(deriv ? c*cm1-s2 : s*cm1);
	__Umm[2][2]=__Umm_coef[2][2]*(deriv ? 2.0*sc : s2);
	__Umm[3][2]=__Umm_coef[3][2]*(deriv ? c*cp1-s2 : s*cp1);
	__Umm[4][2]=__Umm_coef[4][2]*(deriv ? -2.0*cp1*s: std::pow(cp1, 2));
}

void TidalDissipation::calculate_torque_power(const DissipatingBody &body,
		SpinOrbitLockInfo &lock, short body_index,
		Dissipation::Derivative derivative)
{
	for(short lock_dir=-1; lock_dir<=1; ++lock_dir)
		rate_entry(body_index, Dissipation::POWER, derivative, lock_dir)=
			rate_entry(body_index, Dissipation::TORQUEX, derivative,
					lock_dir)=
			rate_entry(body_index, Dissipation::TORQUEZ, derivative,
					lock_dir)=0;
	double spin_frequency=(lock ? lock.spin(__orbital_frequency) :
			body.spin());
	if(body.moment_of_inertia()==0) return;
	for(int m=-2; m<=2; ++m) {
		double m_spin_freq=m*spin_frequency;
		for(int mp=-2; mp<=2; mp+=2) {
			bool locked_term=lock(mp, m);
			PhaseLag::Derivative phase_lag_deriv=PhaseLag::NO_DERIV;
			switch(derivative) {
				case Dissipation::NO_DERIV : case Dissipation::RADIUS :
				case Dissipation::INCLINATION :
					phase_lag_deriv=PhaseLag::NO_DERIV; break;
				case Dissipation::AGE : phase_lag_deriv=PhaseLag::AGE; break;
				case Dissipation::MOMENT_OF_INERTIA : 
				case Dissipation::SPIN_ANGMOM : 
				case Dissipation::SEMIMAJOR :
						   phase_lag_deriv=PhaseLag::FORCING_FREQUENCY;
						   break;
				default :
#ifdef DEBUG
						   assert(false)
#endif
						;
			};
			int mp_ind=mp/2+1, m_ind=m+2;
			double Umm_squared=std::pow(__Umm[m_ind][mp_ind], 2),
				   forcing_freq=(locked_term ? 0 :
						   mp*__orbital_frequency-m_spin_freq);
			for(short lock_dir=(locked_term ? -1 : 0);
					lock_dir<=(locked_term ? 1 : 0);
					lock_dir+=2) {
				lock.lock_direction(lock_dir);
				double mod_phase_lag=body.modified_phase_lag(m, 
						forcing_freq, lock, phase_lag_deriv);
				if(derivative==Dissipation::MOMENT_OF_INERTIA || 
						derivative==Dissipation::SPIN_ANGMOM)
					mod_phase_lag=
						-m*mod_phase_lag
						+
						body.modified_phase_lag(m, forcing_freq, lock,
								PhaseLag::SPIN_FREQUENCY);
				else if(derivative==Dissipation::SEMIMAJOR)
					mod_phase_lag*=mp;
				rate_entry(body_index, Dissipation::TORQUEZ,
						derivative, lock_dir)+=Umm_squared*m*mod_phase_lag;
				rate_entry(body_index, Dissipation::POWER,
						derivative, lock_dir)+=Umm_squared*mp*mod_phase_lag;
				rate_entry(body_index, Dissipation::TORQUEX,
						derivative, lock_dir)+=__Umm[m_ind][mp_ind]*(
							(m>-2 ? 
							 __torque_x_minus_coef[mp_ind]
							 *
							 __Umm[m_ind-1][mp_ind] : 0)+
							(m<2 ? 
							 __torque_x_plus_coef[mp_ind]
							 *
							 __Umm[m_ind+1][mp_ind] : 0));
			}
		}
	}
}

void TidalDissipation::calculate_semimajor_decay(short body_index)
{
	double coef=__mass_product/(2.0*std::pow(__orbital_energy, 2))*\
				AstroConst::G*AstroConst::solar_mass/AstroConst::AU*
				std::pow(AstroConst::day/AstroConst::solar_radius, 2);
	for(short lock_dir=-1; lock_dir<=1; ++lock_dir) {
		for(int fill_deriv=Dissipation::NO_DERIV;
				fill_deriv<Dissipation::NUM_DERIVATIVES; ++fill_deriv)
			rate_entry(body_index, Dissipation::SEMIMAJOR_DECAY, 
					static_cast<Dissipation::Derivative>(fill_deriv),
					lock_dir)=
				coef*rate_entry(body_index, Dissipation::POWER, 
						static_cast<Dissipation::Derivative>(fill_deriv),
						lock_dir);
		rate_entry(body_index, Dissipation::SEMIMAJOR_DECAY,
				Dissipation::SEMIMAJOR, lock_dir)+=rate_entry(body_index,
					Dissipation::POWER, Dissipation::NO_DERIV, lock_dir)/
					__orbital_energy;
		rate_entry(body_index, Dissipation::ORBIT_SPINUP,
				Dissipation::NO_DERIV, lock_dir)=
			1.5*__orbital_frequency/__semimajor*rate_entry(body_index,
					Dissipation::SEMIMAJOR_DECAY, Dissipation::NO_DERIV,
					lock_dir);
	}
}

void TidalDissipation::calculate_inclination_decay(short body_index)
{
	if(
			rate_entry(body_index, Dissipation::TORQUEX,
				Dissipation::NO_DERIV, -1)==0 
			&&
			rate_entry(body_index, Dissipation::TORQUEX,
				Dissipation::NO_DERIV, 0)==0
			&&
			rate_entry(body_index, Dissipation::TORQUEX,
				Dissipation::NO_DERIV, 1)==0
			&&
			rate_entry(body_index, Dissipation::TORQUEZ,
				Dissipation::NO_DERIV, -1)==0 
			&&
			rate_entry(body_index, Dissipation::TORQUEZ,
				Dissipation::NO_DERIV, 0)==0
			&&
			rate_entry(body_index, Dissipation::TORQUEZ,
				Dissipation::NO_DERIV, 1)==0
			&&
			rate_entry(body_index, Dissipation::POWER,
				Dissipation::NO_DERIV, -1)==0 
			&&
			rate_entry(body_index, Dissipation::POWER,
				Dissipation::NO_DERIV, 0)==0
			&&
			rate_entry(body_index, Dissipation::POWER,
				Dissipation::NO_DERIV, 1)==0
			) {
				for(short lock_dir=-1; lock_dir<=1; ++lock_dir) {
					for(int fill_deriv=Dissipation::NO_DERIV;
							fill_deriv<Dissipation::NUM_DERIVATIVES;
							++fill_deriv)
						rate_entry(body_index, Dissipation::INCLINATION_DECAY,
							static_cast<Dissipation::Derivative>(fill_deriv),
							lock_dir)=0;
				}
				return;
			}
	double sin_inclination=std::sin(__inclination[body_index]),
		   cos_inclination=std::cos(__inclination[body_index]),
		   x_coef=1.0/__spin_angular_momentum[body_index]
			   +
			   cos_inclination/__orbital_angular_momentum,
		   z_coef=-sin_inclination/__orbital_angular_momentum;
	for(short lock_dir=-1; lock_dir<=1; ++lock_dir) {
		for(int fill_deriv=Dissipation::NO_DERIV;
				fill_deriv<Dissipation::NUM_DERIVATIVES; ++fill_deriv)
			rate_entry(body_index, Dissipation::INCLINATION_DECAY,
					static_cast<Dissipation::Derivative>(fill_deriv),
					lock_dir)=x_coef*rate_entry(body_index,
						Dissipation::TORQUEX,
						static_cast<Dissipation::Derivative>(fill_deriv),
						lock_dir)
				+
				z_coef*rate_entry(body_index,  Dissipation::TORQUEZ,
						static_cast<Dissipation::Derivative>(fill_deriv),
						lock_dir);
		rate_entry(body_index, Dissipation::INCLINATION_DECAY,
				Dissipation::INCLINATION, lock_dir)-=
			(sin_inclination*rate_entry(body_index, Dissipation::TORQUEX,
										Dissipation::NO_DERIV, lock_dir)
			 +
			 cos_inclination*rate_entry(body_index, Dissipation::TORQUEZ,
				 Dissipation::NO_DERIV, lock_dir)
			)/__orbital_angular_momentum;
		rate_entry(body_index, Dissipation::INCLINATION_DECAY,
				Dissipation::SPIN_ANGMOM, lock_dir)-=rate_entry(body_index,
					Dissipation::TORQUEX, Dissipation::NO_DERIV, lock_dir)/
					std::pow(__spin_angular_momentum[body_index], 2);
		rate_entry(body_index, Dissipation::INCLINATION_DECAY,
				Dissipation::SEMIMAJOR, lock_dir)-=0.5/__semimajor*(
					cos_inclination/__orbital_angular_momentum
					+z_coef);
	}
}

void TidalDissipation::calculate_eccentricity_decay(short)
{
	throw Error::NotImplemented("Eccentric orbits");
}

void TidalDissipation::init(const DissipatingBody &body1, 
		const DissipatingBody &body2, double semimajor, double eccentricity,
		SpinOrbitLockInfo &lock1, SpinOrbitLockInfo &lock2)
{
	__orbital_frequency=orbital_angular_velocity(body1.mass(), body2.mass(),
			semimajor);
	__mass_product=body1.mass()*body2.mass();
	__reduced_mass=__mass_product/(body1.mass()+body2.mass());
	__orbital_energy=orbital_energy(body1.mass(), body2.mass(), semimajor);
	__orbital_angular_momentum=orbital_angular_momentum(body1.mass(),
			body2.mass(), semimajor, eccentricity);
	__semimajor=semimajor;
	__Umm.resize(5, std::valarray<double>(3));
	__dissipation_rate.resize(
			6*Dissipation::NUM_QUANTITIES*Dissipation::NUM_DERIVATIVES, NaN);
	__spin_angular_momentum.resize(2);
	__inclination.resize(2);
#ifdef DEBUG
	__valid=true;
#endif

	__spin_angular_momentum[0]=(lock1 ?
			body1.moment_of_inertia()*lock1.spin(__orbital_frequency) :
			body1.angular_momentum());
	__spin_angular_momentum[1]=(lock2 ? 
			body2.moment_of_inertia()*lock2.spin(__orbital_frequency) :
			body2.angular_momentum());
	__inclination[0]=body1.inclination();
	__inclination[1]=body2.inclination();

	double torque_common=AstroConst::G*AstroConst::solar_mass*
			   			 std::pow(AstroConst::solar_radius/
								  std::pow(semimajor*AstroConst::AU, 2), 3)*
						 AstroConst::Gyr*AstroConst::day,
		   other_mass=body2.mass();
#ifdef DEBUG
	assert(Dissipation::INCLINATION==Dissipation::NUM_DERIVATIVES-1);
#endif
	for(short body_index=0; body_index<2; ++body_index) {
		const DissipatingBody &body=(body_index==0 ? body1 : body2);
		SpinOrbitLockInfo &lock=(body_index==0 ? lock1 : lock2);
		double torque_scale=std::pow(body.radius(),5)*
			std::pow(other_mass, 2)*torque_common,
			spin_frequency=(lock ? lock.spin(__orbital_frequency) : 
					body.spin());
		fill_Umm(body.inclination(), false);
		for(int deriv_ind=Dissipation::NO_DERIV;
				deriv_ind<Dissipation::NUM_DERIVATIVES; ++deriv_ind) {
			Dissipation::Derivative derivative=
				static_cast<Dissipation::Derivative>(deriv_ind);
			switch(derivative) {
				case Dissipation::RADIUS : torque_scale*=5.0/body.radius();
										   break;
				case Dissipation::SEMIMAJOR : torque_scale*=-1.5/semimajor;
											  break;
				case Dissipation::SPIN_ANGMOM :
									torque_scale/=body.moment_of_inertia();
									break;
				case Dissipation::MOMENT_OF_INERTIA :
					torque_scale*=-spin_frequency/body.moment_of_inertia();
					break;
				case Dissipation::INCLINATION :
						fill_Umm(body.inclination(), true);
						break;
				default : ;
			}
			calculate_torque_power(body, lock, body_index,
					derivative);
			for(short lock_dir=-1; lock_dir<=1; ++lock_dir) {
				rate_entry(body_index, Dissipation::POWER, derivative,
						lock_dir)*=torque_scale*__orbital_frequency;
				rate_entry(body_index, Dissipation::TORQUEX, derivative,
						lock_dir)*=torque_scale;
				rate_entry(body_index, Dissipation::TORQUEZ, derivative,
						lock_dir)*=torque_scale;
				if(derivative==Dissipation::SEMIMAJOR) {
					rate_entry(body_index, Dissipation::POWER, derivative,
							lock_dir)-=7.5/semimajor*rate_entry(body_index,
								Dissipation::POWER, Dissipation::NO_DERIV,
								lock_dir);
					rate_entry(body_index, Dissipation::TORQUEX, derivative,
							lock_dir)-=6.0/semimajor*rate_entry(body_index, 
								Dissipation::TORQUEX, Dissipation::NO_DERIV,
								lock_dir);
					rate_entry(body_index, Dissipation::TORQUEZ, derivative,
							lock_dir)-=6.0/semimajor*rate_entry(body_index, 
								Dissipation::TORQUEZ, Dissipation::NO_DERIV,
								lock_dir);
				}
			}
		}
		other_mass=body1.mass();
	}
}

double TidalDissipation::operator()(short body_index, 
		Dissipation::Quantity quantity,
		Dissipation::Derivative derivative,
		short lock_dir) const
{
#ifdef DEBUG
	assert(__valid);
#endif
	if(std::isnan(rate_entry(body_index, quantity, derivative, lock_dir)))
		switch(quantity) {
			case Dissipation::SEMIMAJOR_DECAY :
			case Dissipation::ORBIT_SPINUP : 
				const_cast<TidalDissipation*>(this)->
					calculate_semimajor_decay(body_index);
				break;
			case Dissipation::INCLINATION_DECAY : 
				const_cast<TidalDissipation*>(this)->
					calculate_inclination_decay(body_index);
				break;
			default :
#ifdef DEBUG
				assert(false)
#endif
					;
		}
	return rate_entry(body_index, quantity, derivative, lock_dir);
}
