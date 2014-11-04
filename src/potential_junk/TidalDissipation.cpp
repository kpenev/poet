/**\file
 *
 * \brief Declarations of some methods of the TidalDissipation class.
 *
 * \ingroup StellarSystem_group
 */

#include "TidalDissipation.h"


void TidalDissipation::fill_spin_orbit_harmonics()
{
	__spin_orbit_harmonics.resize(2);
	for(short body_ind=0; body_ind<=1; ++body_ind) {
		__spin_orbit_harmonics[body_ind].resize(4);
		__spin_orbit_harmonics[body_ind][0].set_lock(1, -1, 0);
		__spin_orbit_harmonics[body_ind][0].set_lock(2, -1, 0);
		__spin_orbit_harmonics[body_ind][0].set_lock(2, 1, 0);
		__spin_orbit_harmonics[body_ind][0].set_lock(1, 1, 0);
	}
}

double TidalDissipation::forcing_frequency(short body_index,
		int orbital_multiplier, int spin_multiplier,
		double multiplied_spin)
{
	double nominal=orbital_multiplier*__orbital_frequency-multiplied_spin;
	for(std::vector<SpinOrbitLockInfo>::const_iterator
			harmonic_i=__spin_orbit_harmonics[body_index].begin();
			harmonic_i!=__spin_orbit_harmonics[body_index].end();
			++harmonic_i) {
		if(harmonic_i->term(orbital_multiplier, spin_multiplier)) {
			switch(harmonic_i->lock_direction()) {
				case -1 : return std::min(0.0, nominal);
				case  0 : return 0;
				case  1 : return std::max(0.0, nominal);
				default :
#ifdef DEBUG
						  assert(false)
#endif
							  ;
			}
		}
	}
	std::ostringstream msg;
	msg << "Asked for a forcing frequency for a term (m=" << spin_multiplier
		<< ", mp=" << orbital_multiplier << ") not in the list of tidal "
		"harmonics in TidalDissipation::forcing_frequency!";
	throw Error::BadFunctionArguments(msg.str());
}

void TidalDissipation::calculate_torque_power(const DissipatingBody &body,
		short body_index, Dissipation::Derivative derivative)
{
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
		const SpinOrbitLockInfo &lock1, const SpinOrbitLockInfo &lock2)
{
/*	__orbital_frequency=orbital_angular_velocity(body1.mass(), body2.mass(),
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
	init_harmonics(0, lock1);
	init_harmonics(1, lock2);

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
		const SpinOrbitLockInfo &lock=(body_index==0 ? lock1 : lock2);
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
			calculate_torque_power(body, body_index,
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
#ifdef DEBUG
				if(body_index==0) {
					std::ostringstream msg;
					msg << "NaN result for body: " << body_index << " "
						<< derivative << " derivative of ";
					if(std::isnan(rate_entry(body_index, Dissipation::POWER,
									derivative, lock_dir)))
						throw Error::Runtime(msg.str()+"POWER");
					if(std::isnan(rate_entry(body_index,Dissipation::TORQUEX,
									derivative, lock_dir)))
						throw Error::Runtime(msg.str()+"TORQUEX");
					if(std::isnan(rate_entry(body_index,Dissipation::TORQUEZ,
									derivative, lock_dir)))
						throw Error::Runtime(msg.str()+"TORQUEZ");
				}
#endif
			}
		}
		other_mass=body1.mass();
	}*/
}

double TidalDissipation::operator()(short body_index, 
		Dissipation::Quantity quantity,
		Dissipation::Derivative derivative,
		short lock_dir) const
{
#ifdef DEBUG
	assert(__valid);
#endif
	if(std::isnan(rate_entry(body_index, quantity, derivative, lock_dir))
		&& quantity!=Dissipation::POWER && quantity!=Dissipation::TORQUEX &&
			quantity!=Dissipation::TORQUEZ)
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

double TidalDissipation::operator()(short body_index, 
		Dissipation::Quantity quantity, double above_fraction,
		Dissipation::Derivative derivative) const
{
	return (*this)(body_index, quantity, derivative, 0) +
		above_fraction*(*this)(body_index, quantity, derivative, 1) + 
		(1.0-above_fraction)*(*this)(body_index, quantity, derivative, -1);
}

void TidalDissipation::init_harmonics(short body_index,
		const SpinOrbitLockInfo &lock)
{
	int orbital_frequency_multiplier=lock.orbital_frequency_multiplier(),
		spin_frequency_multiplier=lock.spin_frequency_multiplier();
	if(orbital_frequency_multiplier==-spin_frequency_multiplier)
		__spin_orbit_harmonics[body_index][0].lock_direction(
				lock.lock_direction());
	else __spin_orbit_harmonics[0][body_index].lock_direction(
			(spin_frequency_multiplier*
			 (orbital_frequency_multiplier+spin_frequency_multiplier)>0
			 ? 1 : -1));
	if(2*spin_frequency_multiplier==-orbital_frequency_multiplier)
		__spin_orbit_harmonics[body_index][1].lock_direction(
				lock.lock_direction());
	else __spin_orbit_harmonics[body_index][1].lock_direction( 
			(spin_frequency_multiplier*
			 (2*spin_frequency_multiplier+orbital_frequency_multiplier)>0
			 ? 1 : -1));
	if(2*spin_frequency_multiplier==orbital_frequency_multiplier)
		__spin_orbit_harmonics[body_index][2].lock_direction(
				lock.lock_direction());
	else __spin_orbit_harmonics[body_index][2].lock_direction(
			(spin_frequency_multiplier*
			 (2*spin_frequency_multiplier-orbital_frequency_multiplier)>0
			 ? 1 : -1));
	if(spin_frequency_multiplier==orbital_frequency_multiplier)
		__spin_orbit_harmonics[body_index][3].lock_direction(
				lock.lock_direction());
	else __spin_orbit_harmonics[body_index][3].lock_direction(
			(spin_frequency_multiplier*
			 (spin_frequency_multiplier-orbital_frequency_multiplier)>0
			 ? 1 : -1));
}

void TidalDissipation::init_harmonics(short body_index,
		double spin_to_orbital_ratio)
{
	for(unsigned i=0; i<__spin_orbit_harmonics.size(); ++i) {
#ifdef DEBUG
		assert(__spin_orbit_harmonics[body_index][i].
				orbital_frequency_multiplier()
				!=
				spin_to_orbital_ratio
				*
				__spin_orbit_harmonics[i][body_index].
				spin_frequency_multiplier());
#endif
		if(__spin_orbit_harmonics[body_index][i].
				orbital_frequency_multiplier()
				>
				spin_to_orbital_ratio
				*
				__spin_orbit_harmonics[body_index][i].
				spin_frequency_multiplier())
			__spin_orbit_harmonics[body_index][i].lock_direction(1);
		else __spin_orbit_harmonics[body_index][i].lock_direction(-1);
	}
}