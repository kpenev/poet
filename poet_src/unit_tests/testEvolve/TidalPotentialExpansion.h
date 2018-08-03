/**\file
 *
 * \brief Declare an interface for evaluating the expansion of the tidal
 * potential.
 *
 * \ingroup UnitTests_group
 */

#ifndef __UNIT_TESTS_TIDAL_POTENTIAL_EXPANSION_H
#define __UNIT_TESTS_TIDAL_POTENTIAL_EXPANSION_H

#include "EccentricOrbit.h"

#include "../../Evolve/TidalPotentialTerms.h"

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cassert>

namespace Evolve {
    ///Evaluate the tidal potential using the expansion.
    class TidalPotentialExpansion {
    private:
        ///\brief The coefficients of the expansion of the tidal potential.
        ///(\f$ \mathcal{U}_{m,m'} \f$)
        TidalPotentialTerms __expansion_coef;

        double
            ///The mass of the tidally perturbed object in solar masses.
            __primary_mass,

            ///The mass of the perturber object in solar masses.
            __secondary_mass,

            ///The semimajor axis of the orbit in solar radii.
            __semimajor,
            
            ///The eccentricity of the orbit
            __eccentricity,

            ///The angle between the orbital angular momentum and the spin
            ///angular momentum of the primary.;
            __inclination,

            ///90 degrees less than the angle from
            /// \f$ \hat{y} = \hat{S} \times \hat{L} \f$
            ///to the direction of periapsis in radians.
            __arg_of_periapsis;

        ///\brief Return a single tidal term:
        /// \f$ \Sum_{m=-2}^{2} \mathcal{U}_{m,m'} \rho'%2 Y_{2,m}(\theta', \phi')\exp(-im'\Omega t) \f$
        double tidal_term(
            ///The value of m' in the expression evaluated by this function.
            int mprime,
            
            ///See same name argument to evaluate_spherical_coords()
            double radial_distance,

            ///See same name argument to evaluate_spherical_coords()
            double azimuthal_angle,

            ///See same name argument to evaluate_spherical_coords()
            double polar_angle,

            ///The orbital phase of the secondary body
            ///(\f$ 2 \pi n \f$ is periapsis).
            double orbital_phase
        );
    public:
        TidalPotentialExpansion(
            ///See same name argument to EccentricOrbit.
            double primary_mass=Core::NaN,

            ///See same name argument to EccentricOrbit.
            double secondary_mass=Core::NaN,

            ///See same name argument to EccentricOrbit.
            double semimajor=Core::NaN,

            ///See same name argument to EccentricOrbit.
            double eccentricity=Core::NaN,

            ///See same name argument to TidalPotential.
            double inclination=Core::NaN,

            ///See same name argument to TidalPotential.
            double arg_of_periapsis=Core::NaN
        ) :
            __primary_mass(primary_mass),
            __secondary_mass(secondary_mass),
            __semimajor(semimajor),
            __eccentricity(eccentricity),
            __inclination(inclination),
            __arg_of_periapsis(arg_of_periapsis)
        {}

        ///Return the tidal potential at a specific position in polar
        ///coordinates and time in seconds.
        double evaluate_spherical_coords(
            ///The radial distance from the origin of the point to where to
            ///evaluate the tidal potential.
            double radial_distance,

            ///The azimuthal angle of the point to where to evaluate the tidal
            ///potential, should be in the range of \f$ [0, 2\pi) \f$.
            double azimuthal_angle,

            ///The polar angle of the point to where to evaluate the tidal
            ///potential, should be in the range of \f$ [0, \pi] \f$.
            double polar_angle,

            ///The time in second since periastron passage. It is perfectly
            ///valid to pass values bigger than the orbital period.
            double time
        );

        ///Return the tidal potential at a specific position and time in SI.
        template<class POSITION_TYPE>
            double operator()(
                ///The position to evaluate the potential at in a coordinate
                ///system centered on the primary body with
                /// \f$ \hat{z} = \hat{S} \f$,
                /// \f$ \hat{y} = \hat{S} \times \hat{L} \f$.
                ///
                ///Must provide indexing with indices 0, 1, 2 for the three
                ///components.
                const POSITION_TYPE &position,

                ///The time in days when to evalutae the tidal potential. The
                ///system is in periapsis at time = 0.
                double time
            );

        void set_eccentricity_order(unsigned e_order)
        {__expansion_coef.change_e_order(e_order);}
    }; //End TidalPotentialExpansion class.

    template<class POSITION_TYPE>
        double TidalPotentialExpansion::operator()(
            const POSITION_TYPE &position,
            double time
        )
        {
            double radial_distance = position.norm(),
                   polar_angle = (
                       radial_distance == 0
                       ? 0
                       : std::acos(position[2]/radial_distance)
                   ),
                   azimuthal_angle = std::atan2(position[1], position[0]);

            if(azimuthal_angle < 0)
                azimuthal_angle += 2.0 * M_PI;

            return evaluate_spherical_coords(
                radial_distance,
                azimuthal_angle,
                polar_angle,
                time
            );
        }
} //End Evolve namespace.
#endif
