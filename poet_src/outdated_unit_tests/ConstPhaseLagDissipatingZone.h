/**\file
 *
 * \brief Declares a dissipating zone with constant but different phase lags
 * for each tidal component.
 *
 * \ingroup UnitTests_group
 */

#ifndef __CONST_PHASE_LAG_DISSIPATING_ZONE_H
#define __CONST_PHASE_LAG_DISSIPATING_ZONE_H

#include "../BinarySystem/DissipatingZone.h"
#include <map>
#include <utility>

inline bool operator<(const std::pair<int, int> &v1, 
					  const std::pair<int, int> &v2)
{
	return v1.first<v2.first || (v1.first==v2.first && v1.second<v2.second);
}

class Lags : public std::map< std::pair<int, int>, double> {
public:
	///A reference to the requested element (creates it if it does not
	///exist).
	double &operator()(int m, int mp);

	///A copy of the requested element (error if it does not exist).
	double operator()(int m, int mp) const;
};

class TestingDissipatingZone : public DissipatingZone {
public:
	///Should place enough informatino to identify the zone to the given
	///stream.
	virtual void describe(std::ostream &os) const=0;
};

class ConstPhaseLagDissipatingZone : public TestingDissipatingZone {
private:
	Lags __lags;

	std::vector<double> __moment_of_inertia, __radius, __mass;

	bool __flip_sign;
public:
	ConstPhaseLagDissipatingZone(const Lags &lags, double moment_of_inertia,
			double radius, double mass, bool flip_sign=false,
			double moment_of_inertia_deriv=0,
			double radius_deriv=0, double mass_deriv=0)
		: __lags(lags), __moment_of_inertia(3), __radius(3), __mass(3),
		__flip_sign(flip_sign)
	{
		__moment_of_inertia[0]=moment_of_inertia; 
		__moment_of_inertia[1]=moment_of_inertia_deriv;
		__moment_of_inertia[2]=0;
		__radius[0]=radius;
		__radius[1]=radius_deriv;
		__radius[2]=0;
		__mass[0]=mass;
		__mass[1]=mass_deriv;
		__mass[2]=0;
	}

	///\brief Should return true iff the given term is presently locked.
	virtual bool locked(int orbital_frequency_multiplier,
			int spin_frequency_multiplier) const {return false;}

	///\brief Should return the tidal phase lag time the love number for the
	///given tidal term (or one of its derivatives).
	///
	///In case the specified term is in a lock, it should return the phase
	///lag for the case of the spin frequency approaching the lock from
	///below. The lag for spin frequency approaching from above should be
	///written to above_lock_value. If the term is not locked 
	///leave above_lock_value untouched.
	virtual double modified_phase_lag(
			///The multiplier of the orbital frequency in the
			///expression for the forcing frequency.
			int orbital_frequency_multiplier,

			///The multiplier of the spin frequency in the
			///expression for the forcing frequency.
			int spin_frequency_multiplier,
			
			///The current forcinc frequency in rad/day
			double forcing_frequency,

			///The return value should be either the phase lag itself
			///(NO_DERIV) or its derivative w.r.t. the specified quantity.
			Dissipation::Derivative deriv,

			///If the lag of a locked term is calculated this should be set
			///to the lag assuming the spin frequency is just above the lock.
			///Otherwise, leave untouched.
			double &above_lock_value) const
	{
		if(deriv!=Dissipation::NO_DERIV) return 0;
		if(lock_held()(orbital_frequency_multiplier,
					   spin_frequency_multiplier)) {
			above_lock_value=(spin_frequency_multiplier>0 ? -1 : 1)
							 *__lags(spin_frequency_multiplier,
									 orbital_frequency_multiplier);
			return -above_lock_value;
		}
		return (__flip_sign && forcing_frequency<0 ? -1 : 1)*
			   __lags(spin_frequency_multiplier,
					  orbital_frequency_multiplier);
	}

	Lags &lags() {return __lags;}

	double moment_of_inertia(int deriv) const
	{return __moment_of_inertia.at(deriv);}

	double moment_of_inertia(double, int deriv) const
	{return moment_of_inertia(deriv);}

	double love_coefficient(int, int, Dissipation::Derivative) const
	{return 0;}

	double outer_radius(int deriv) const
	{
		return __radius.at(deriv);
	}

	double outer_radius(double age, int deriv) const
	{
		return outer_radius(deriv);
	}

	double outer_mass(int deriv) const 
	{
		return __mass.at(deriv);
	}

	double outer_mass(double age, int deriv) const
	{
		return outer_mass(deriv);
	}

	bool dissipative() const {return true;}

	void describe(std::ostream &os) const;

	virtual ~ConstPhaseLagDissipatingZone() {}
};


#endif
