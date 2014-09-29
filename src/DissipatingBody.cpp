void DissipatingBody::angular_momentum_transfer(
		const DissipatingZone &outer_zone,
		const DissipatingZone &inner_zone,
		Eigen::Vector3D &outer_angmom_gain,
		Eigen::Vector3D &inner_angmom_gain,
		Dissipation::Derivative deriv, bool with_respect_to_outer)
{
#ifdef DEBUG
	assert(deriv==Dissipation::NO_DERIV || deriv==Dissipation::INCLINATION
			|| deriv==Dissipation::PERIAPSIS);
#endif
	double dm_dt=inner_zone.outer_mass(1),
		   lost_spin=(dm_dt>0 ? outer_zone : inner_zone).spin_frequency(),
		   angmom_transfer=-2.0/3.0*std::pow(inner_zone.outer_radius(), 2)*
			   lost_spin*dm_dt;
	if(dm_dt>0) {
		outer_angmom_gain[0]=outer_angmom_gain[1]=0;
		outer_angmom_gain[2]=
			(deriv==Dissipation::NO_DERIV ? angmom_transfer : 0);
		inner_angmom_gain=-zone_to_zone_transform(outer_zone, inner_zone,
				Eigen::Vector3d(0, 0, angmom_transfer), deriv, 
				with_respect_to_outer);
	} else {
		inner_angmom_gain[0]=inner_angmom_gain[1]=0;
		inner_angmom_gain[2]=
			(deriv==Dissipation::NO_DERIV ? angmom_transfer : 0);
		inner_angmom_gain=-zone_to_zone_transform(inner_zone, outer_zone,
				Eigen::Vector3d(0, 0, angmom_transfer), deriv, 
				!with_respect_to_outer);
	}
}

Eigen::Vector3d DissipatingBody::angular_momentum_transfer_from_top(
		unsigned zone_index, Dissipation::Derivative deriv, 
		bool with_respect_to_outer)
{
#ifdef DEBUG
	assert(zone_index>0);
#endif
	DissipatingZone &this_zone=zone(zone_index),
					&zone_above=zone(zone_index-1);
	double scaling, dm_dt=this_zone.outer_mass(1);
	if(deriv==Dissipation::AGE) {
		scaling=2.0*this_zone.outer_radius(1)/this_zone.outer_radius(0)
			    + this_zone.outer_mass(2)/dm_dt;
	} else if(deriv==Dissipation::SPIN_FREQUENCY || 
			  deriv==Dissipation::MOMENT_OF_INERTIA ||
			  deriv==Dissipation::SPIN_ANGMOM) {
		if((with_respect_to_outer && dm_dt<=0) ||
				(!with_respect_to_outer && dm_dt>=0)) 
			return Eigen::Vector3d(0, 0, 0);
		else if(deriv==Dissipation::SPIN_FREQUENCY)
			scaling=1.0/(with_respect_to_outer ? zone_above.spin_frequency()
									           : this_zone.spin_frequency());
		else if(deriv==Dissipation::MOMENT_OF_INERTIA)
			scaling=-1.0/(with_respect_to_outer
						  ? zone_above.moment_of_inertia()
						  : this_zone.moment_of_inertia());
		else scaling=1.0/(with_respect_to_outer
						  ? zone_above.angular_momentum()
						  : this_zone.angular_momentum());
	} else if(deriv==Dissipation::INCLINATION || 
			  deriv==Dissipation::PERIAPSIS) {
		if(dm_dt<=0) return Eigen::Vector3d(0, 0, 0);
		else {
			Eigen::Vector3d dummy, result;
			angular_momentum_transfer(zone_above, this_zone, dummy, result,
					deriv, with_respect_to_outer);
			return result;
		}
	}
#ifdef DEBUG
	else assert(false);
#endif
	return scaling*__angular_momentum_transfer[zone_index-1][1];
}

Eigen::Vector3d DissipatingBody::angular_momentum_transfer_from_bottom(
		unsigned zone_index, Dissipation::Derivative deriv, 
		bool with_respect_to_inner)
{
#ifdef DEBUG
	assert(zone_index<number_zones()-1);
#endif
	DissipatingZone &this_zone=zone(zone_index),
					&zone_below=zone(zone_index+1);
	double scaling, dm_dt=zone_below.outer_mass(1);
	if(deriv==Dissipation::AGE) {
		scaling=2.0*zone_below.outer_radius(1)/zone_below.outer_radius(0)
			    + zone_below.outer_mass(2)/dm_dt;
	} else if(deriv==Dissipation::SPIN_FREQUENCY || 
			  deriv==Dissipation::MOMENT_OF_INERTIA ||
			  deriv==Dissipation::SPIN_ANGMOM) {
		if((with_respect_to_inner && dm_dt>=0) ||
				(!with_respect_to_inner && dm_dt<=0))
			return Eigen::Vector3d(0, 0, 0);
		else if(deriv==Dissipation::SPIN_FREQUENCY)
			scaling=1.0/(with_respect_to_inner ? zone_below.spin_frequency()
					                           : this_zone.spin_frequency());
		else if(deriv==Dissipation::MOMENT_OF_INERTIA)
			scaling=-1.0/(with_respect_to_inner
					      ? zone_below.moment_of_inertia()
						  : this_zone.moment_of_inertia());
		else scaling=1.0/(with_respect_to_inner
				          ? zone_below.angular_momentum()
						  : this_zone.angular_momentum());
	} else if(deriv==Dissipation::INCLINATION ||
			  deriv==Dissipation::PERIAPSIS) {
		if(dm_dt>=0) return Eigen::Vector3d(0, 0, 0);
		else {
			Eigen::Vector3d dummy, result;
			angular_momentum_transfer(this_zone, zone_below, result, dummy,
					deriv, !with_respect_to_inner);
		}
	}
#ifdef DEBUG
	else assert(false)
#endif
	return scaling*__angular_momentum_transfer[zone_index][0];
}

Eigen::Vector3d DissipatingBody::angular_momentum_transfer_to_zone(
			unsigned zone_index, Dissipation::Derivative deriv,
			int deriv_zone)
{
	if(deriv==Dissipation::ORBITAL_FREQUENCY
			|| deriv==Dissipation::ECCENTRICITY
			|| deriv==Dissipation::RADIUS
			|| deriv==Dissipation::SEMIMAJOR) return Eigen::Vector3d(0,0,0);
	Eigen::Vector3d result(0, 0, 0);
	if(zone_index>0 && (deriv==Dissipation::NO_DERIV || 
			deriv==Dissipation::AGE || deriv_zone<=0)) 
		result=angular_momentum_transfer_from_top(zone_index, deriv, 
				deriv_zone<0);
	if(zone_index<num_zones()-1 && (deriv==Dissipation::NO_DERIV || 
			deriv==Dissipation::AGE || deriv_zone>=0))
		result+=angular_momentum_transfer_from_bottom(zone_index, deriv, 
				deriv_zone>0);
	return result;
}

void DissipatingBody::set_orbit(double companion_mass,
		double semimajor, double eccentricity)
{
	double torque_norm=std::pow(companion_mass/std::pow(semimajor, 3), 2)
					   *std::pow(radius(), 5)
					   *AstroConst::G*AstroConst::day*AstroConst::Gyr
					   *AstroConst::solar_mass
					   /std::pow(AstroConst::solar_radius, 3),
		   orbital_frequency=orbital_angular_velocity(companion_mass, mass(),
				   semimajor, false)
		   power_norm=tidal_torque_norm*orbital_frequency,
		   dangular_velocity_da=-1.5*orbital_frequency/semimajor;
				
	__tidal_torques.resize(number_zones());
	__tidal_power=0;
	__angular_momentum_transfer(number_zones()-1);
	for(zone_index=0; zone_index<number_zones(); ++zone_index) {
		DissipatingZone &current_zone=zone(zone_index);
		current_zone.set_orbit(orbital_frequency, eccentricity);
		if(zone_index<number_zones()-1) {
			__angular_momentum_transfer[zone_index].resize(2);
			angular_momentum_transfer(zone_index, zone_index+1, 
					__angular_momentum_transfer[zone_index][0],
					__angular_momentum_transfer[zone_index][1]);
		}
		for(above=false; !above; above=true) {
			std::valarray<Eigen::Vector3d> &tidal_torque=
				(above ? __tidal_torques_above : __tidal_torques_below)
				[zone_index];
			tidal_torque.resize(Dissipation::NUM_DERIVATIVES)
			tidal_power=(above ? __tidal_power_above : __tidal_power_below);
			for(Dissipation::Derivative deriv=Dissipation::NO_DERIV;
				deriv<Dissipation::END_DIMENSIONLESS_DERIV; ++deriv) {
				tidal_power[deriv]+=current_zone.tidal_power(deriv, true);
				tidal_torque[deriv][0]=
					current_zone.tidal_torque_x(deriv, above);
				tidal_torque[deriv][1]=
					current_zone.tidal_torque_y(deriv, above);
				tidal_torque[deriv][2]=
					current_zone.tidal_torque_z(deriv, above);
			}
			tidal_torque[Dissipation::ORBITAL_FREQUENCY]+=
				4.0/orbital_frequency*tidal_torque[Dissipation::NO_DERIV];
			tidal_torque[Dissipation::RADIUS]=
				5.0/radius()*tidal_torque[Dissipation::NO_DERIV];
			tidal_power[Dissipation::ORBITAL_FREQUENCY]+=
				5.0/orbital_frequency*tidal_power[Dissipation::NO_DERIV];
			tidal_power[Dissipation::RADIUS]+=
				5.0/radius()*tidal_power[Dissipation::NO_DERIV];
			tidal_torque[Dissipation::MOMENT_OF_INERTIA]=
				-current_zone.spin_frequency()
				/current_zone.moment_of_inertia()
				*tidal_torque[Dissipation::SPIN_FREQUENCY];
			tidal_torque[Dissipation::SPIN_ANGMOM]=
				tidal_torque[Dissipation::SPIN_FREQUENCY]
				/current_zone.moment_of_inertia();
			tidal_torque[Dissipation::SEMIMAJOR]=
				dangular_velocity_da
				*tidal_torque[Dissipation::ORBITAL_FREQUENCY];
		}
	}
	__tidal_power_above*=power_norm;
	__tidal_power_below*=power_norm;
	__tidal_torques_above*=torque_norm;
	__tidal_torques_below*=torque_norm;
}

const Eigen::Vector3D &external_torque(unsigned zone_index,
		Dissipation::Derivative deriv, int deriv_zone)
{
#ifdef DEBUG
	assert(zone_index<number_zones());
	assert(static_cast<int>(zone_index)+deriv_zone>0);
	assert(static_cast<int>(zone_index)+deriv_zone<number_zones());
#endif
	DissipatingZone &this_zone=zone(zone_index);
	Eigen::Vector3D result(0, 0, 0);
	if(zone_index==0 && deriv_zone==0)
		result[2]=angular_momentum_loss(deriv);
	result+=angular_momentum_transfer_to_zone(zone_index, deriv, deriv_zone);
	if(zone_index<number_zones()-1 && 
			(!zone_specific(deriv) || deriv_zone>=0))
		result+=angular_momentum_coupling(zone_index, deriv, deriv_zone==0);
	if(zone_index>0 && (!zone_specific(deriv) || deriv_zone<=0)) {
		result+=zone_to_zone_transform(
				zone(zone_index-1), this_zone,
				angular_momentum_coupling(zone_index-1, deriv,
					                      deriv_zone<0));
		if(deriv==Dissipation::INCLINATION || deriv==Dissipation::PERIAPSIS)
			result+=zone_to_zone_transform(
					zone(zone_index-1), this_zone,
					angular_momentum_coupling(zone_index-1), deriv, 
					deriv_zone<0);
	}
	return result;
}
