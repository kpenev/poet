/**\file
 *
 * \brief Declares a class for the single zone of LockedPlanet objects.
 *
 * \ingroup Planet_group
 */

#ifndef __LOCKED_PLANET_ZONE_H
#define __LOCKED_PLANET_ZONE_H

#include "../Core/SharedLibraryExportMacros.h"
#include "../Evolve/DissipatingZone.h"
#include "../Core/AstronomicalConstants.h"
#include "../Evolve/BrokenPowerlawPhaseLagZone.h"

namespace Planet {

    ///\brief The only zone of a LockedPlanet.
    ///
    ///\ingroup Planet_group
    class LIB_LOCAL PlanetZone :
        virtual public Evolve::BrokenPowerlawPhaseLagZone {
    private:
        double
            ///See mass argument to constructor
            __mass,

            ///See radius argument to constructor
            __radius,

            ///See moment_of_inertia_factor to constructor.
            __moment_of_inertia_factor;

    public:
        PlanetZone(
            ///The mass of the planet in solar masses.
            double mass,

            ///The radius of the planet is solar radii.
            double radius,

            ///The factor to multiply mass * radius^2 to get the moment of
            ///inertia.
            double moment_of_inertia_factor
        ) :
            __mass(mass),
            __radius(radius),
            __moment_of_inertia_factor(moment_of_inertia_factor)
        {}

        ///See DissipatingZone::love_coefficient(), always zero.
        double love_coefficient(int,
                                int,
                                Evolve::Dissipation::QuantityEntry) const
        {return 0;}

        ///Tiny value (\f$10^{-6} M R^2\f$).
        double moment_of_inertia(
            ///What to return:
            /// - 0 The moment of inertia in \f$M_\odot R_\odot^2\f$
            /// - 1 The rate of change of the moment of inertia in
            ///     \f$M_\odot R_\odot^2/Gyr\f$
            /// - 2 The second derivative in \f$M_\odot R_\odot^2/Gyr^2\f$
            int deriv_order=0
        ) const
        {
            return (
                deriv_order==0
                ? __moment_of_inertia_factor * __mass * std::pow(__radius, 2)
                : 0
            );
        }

        double moment_of_inertia(double, int deriv_order=0) const
        {return moment_of_inertia(deriv_order);}

        ///The radius of the planet.
        double outer_radius(int deriv_order=0) const
        {return (deriv_order==0 ? __radius : 0);}

        ///Same as outer_radius(int) but accept age argument (ignored).
        double outer_radius(double, int deriv_order=0) const
        {return (deriv_order==0 ? __radius : 0);}

        ///The mass of the planet.
        double outer_mass(int deriv_order=0) const
        {return (deriv_order==0 ? __mass : 0);}

        ///Same as outer_mass(int) but accept age argument (ignored).
        double outer_mass(double, int deriv_order=0) const
        {return (deriv_order==0 ? __mass : 0);}

        ///\brief Calls the usual DissipatingZone::configure but with zero
        ///inclination and periapsis.
        void configure(
            ///Is this the first time configure() is invoked?
            bool initialize,

            ///The age to set the zone to.
            double age,

            ///The angular velocity of the orbit in rad/day.
            double orbital_frequency,

            ///The eccentricity of the orbit
            double eccentricity,

            ///The absolute value of the angular momentum of the orbit.
            double orbital_angmom,

            ///The angular momentum/frequency of the spin of the zone
            ///(ignored).
            double spin,

            ///The inclination of the zone relative to the orbit.
            double inclination,

            ///The argument of periapsis of the orbit in the equatorial
            ///planet of the zone.
            double periapsis,

            ///Is spin_angmom angular momentum of freuqency? (ignored)
            bool spin_is_frequency,

            ///If specified, only a single tidal term is included in the tidal
            ///rates. The two values specify the spin and orbital frequency
            ///multipliers of the term to add (in that order).
            std::pair<int, int> *single_term=NULL
        )
        {
            BrokenPowerlawPhaseLagZone::configure(
                initialize,
                age,
                orbital_frequency,
                eccentricity,
                orbital_angmom,
                (dissipative() ? spin : orbital_frequency),
                (dissipative() ? inclination : 0),
                (dissipative() ? periapsis : 0),
                (dissipative() ? spin_is_frequency : true),
                single_term
            );
        }

    }; //End LockedPlanetZone class.

}//End Planet namespace.

#endif
