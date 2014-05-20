/**\file 
 *
 * \brief Declares the DissipatingBody class.
 *
 * \ingroup StellarSystem_group
 */

///\brief A base class for any body contributing to tidal dissipation.
///
///\ingroup StellarSystem_group
class DissipatingBody {
private:
	///The spin frequency of the body in rad/day.
	double __spin_frequency,

		   ///Angle between the spin and orbital angular momentum in radians.
		   __inclination;
public:
	///\brief A function defining the dissipation efficiency of the body.
	///
	///In our formalism this function should return
	/// \f$\Delta_{m,m'}'\equiv\kappa_{m,m'}\sin(\Delta_{m,m'})\f$.
	///For small dissipation efficiency, this is the phase lag times the love
	///number. It must satisfy \f$\Delta_{m,m'}=\Delta_{-m,-m'}\f$.
	virtual double modified_phase_lag(
			///The \f$m\f$ index.
			int m, 

			///The forcing frequency (\f$m'\Omega-m\Omega_s\f$ in Lai, 2012)
			///in rad/day.
			double forcing_frequency,
			
			///If this function is discontinuous at zero and an exactly zero
			///forcing frequency is encountered, this flag determines if the
			///zero should be interpreted as an infinitesimal positive or
			///negative amount above. It is safe to ignore this flag if the
			///function is continuous at zero forcing frequency.
			short forcing_sign) const=0;

	///The radius of this body in \f$R_\odot\f$.
	virtual double radius() const =0;

	///The mass of this body in \f$M_\odot\f$.
	virtual double mass() const =0;

	///The spin frequency of the body in rad/day.
	double spin_frequency() const {return __spin_frequency};

	///Reference to the spin frequency of the body in rad/day.
	double &spin_frequency() {return __spin_frequency};

	///Angle between the spin and orbital angular momentum in radians.
	double inclination() const {return __inclination};

	///\brief Reference to the angle between the spin and orbital angular 
	///momentum in radians.
	double &inclination() {return __inclination};
};
