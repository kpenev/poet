/**\file
 *
 * \brief Declare an interface for working with eccentric orbits.
 *
 * \ingroup UnitTests_group
 */

#ifndef __UNIT_TESTS_ECCENTRIC_ORBIT_H
#define __UNIT_TESTS_ECCENTRIC_ORBIT_H

#include "../../Core/Common.h"
#include "../../Core/AstronomicalConstants.h"

#include "Eigen/Dense"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include <iostream>
#include <iomanip>

namespace Evolve {
    ///Basic description of two bodies in an eccentric orbit.
    class EccentricOrbit {
    private:
        double
            ///The mass of the tidally perturbed object in solar masses.
            __primary_mass,

            ///The mass of the perturber object in solar masses.
            __secondary_mass,

            ///The semimajor axis of the orbit in solar radii.
            __semimajor,

            ///The eccentricity of the orbit
            __eccentricity;

        ///The reduced mass of the system in solar masses.
        double reduced_mass() const;

    public:
        ///Create an eccentric orbit.
        EccentricOrbit(
            ///See __primary_mass attribute.
            double primary_mass=Core::NaN,

            ///See __secondary_mass attribute.
            double secondary_mass=Core::NaN,

            ///See __semimajor attribute.
            double semimajor=Core::NaN,

            ///See __eccentricity attribute.
            double eccentricity=Core::NaN
        ) :
            __primary_mass(primary_mass),
            __secondary_mass(secondary_mass),
            __semimajor(semimajor),
            __eccentricity(eccentricity)
        {}

        ///The semimajor axis of the system.
        double primary_mass() const {return __primary_mass;}

        ///A mutable reference to the  semimajor axis of the system.
        double &primary_mass() {return __primary_mass;}

        ///The semimajor axis of the system.
        double secondary_mass() const {return __secondary_mass;}

        ///A mutable reference to the  semimajor axis of the system.
        double &secondary_mass() {return __secondary_mass;}

        ///The semimajor axis of the system.
        double semimajor() const {return __semimajor;}

        ///A mutable reference to the  semimajor axis of the system.
        double &semimajor() {return __semimajor;}

        ///The semimajor axis of the system.
        double eccentricity() const {return __eccentricity;}

        ///A mutable reference to the  semimajor axis of the system.
        double &eccentricity() {return __eccentricity;}

        ///The eccentric anomaly that corresponds to the given orbital phase
        ///angle.
        double eccentric_anomaly(
            ///The phase angle of the orbit (\f$ \Omega t \f$). Zero is at
            ///periapsis
            double orbital_phase
        ) const;

        ///\brief Secondary position vector in a coordinate system centered on
        ///the primary, with \f$ \hat{z} = \hat{L} \f$ and
        /// \f$ \hat{y} = \hat{S} \times \hat{L} \f$
        Eigen::Matrix<long double, 3, 1> secondary_position(
            ///See same name argument to eccentric_anomaly().
            double orbital_phase
        ) const;

        ///\brief The magnitude of the orbital angular momentum of the system in
        /// \f$ M_\odot R_\odot^2 \mathrm{rad} / \mathrm{Gyr} \f$
        double orbital_angmom() const;

        ///\brief The magnitude of the orbital energy of the system in
        /// \f$ M_\odot R_\odot^2 / \mathrm{Gyr}^2 \f$
        double orbital_energy() const;

        ///\brief The orbital period of the system in days.
        double orbital_period() const;
    }; //End EccentricOrbit class
}//End Evolve namespace

#endif
