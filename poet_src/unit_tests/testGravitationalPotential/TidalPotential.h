/**\file
 *
 * \brief Declare an interface for calculating the tidal potential due to a
 * companion in an eccentric orbit.
 *
 * \ingroup UnitTests_group
 */

#include "EccentricOrbit.h"
#include <cmath>

namespace testGravitationalPotential {
    ///Calculate the tidal potential over one component of an eccentric binary.
    class TidalPotential {
    private:
        ///The binary orbit.
        EccentricOrbit __orbit;

        double
            ///The angle from \f$ \hat{y} = \hat{S} \times \hat{L} \f$ to the
            ///direction of periapsis in radians.
            __arg_of_periapsis;
    public:
        ///\brief Define the boundary for which to calculate the tidal
        ///potential.
        TidalPotential(
            ///See same name argument to EccentricOrbit.
            double primary_mass=Core::NaN,

            ///See same name argument to EccentricOrbit.
            double secondary_mass=Core::NaN,

            ///See same name argument to EccentricOrbit.
            double semimajor=Core::NaN,

            ///See same name argument to EccentricOrbit.
            double eccentricity=Core::NaN,
            
            ///See __arg_of_periapsis attribute.
            double arg_of_periapsis=Core::NaN
        ) :
            __orbit(primary_mass, secondary_mass, semimajor, eccentricity),
            __arg_of_periapsis(arg_of_periapsis)
        {}

        ///An unmutable reference to the binary orbit.
        const EccentricOrbit &orbit() const {return __orbit;}

        ///Mutable reference to the binary orbit.
        EccentricOrbit &orbit() {return __orbit;}

        ///Return the tidal potential at a specific position and time in SI.
        template<class POSITION_TYPE>
            double operator()(
                ///The position to evaluate the potential at in a coordinate
                ///system centered on the primary body with
                /// \f$ \hat{z} = \hat{L} \f$,
                /// \f$ \hat{y} = \hat{S} \times \hat{L} \f$.
                ///
                ///Must provide indexing with indices 0, 1, 2 for the three
                ///components.
                const POSITION_TYPE &position,

                ///The time in days when to evalutae the tidal potential. The
                ///system is in periapsis at time = 0.
                double time
            ) const;
    };

    template<class POSITION_TYPE>
        double TidalPotential::operator()(const POSITION_TYPE &position,
                                          double time) const
        {
            Eigen::Vector3d secondary_position = __orbit.secondary_position(
                2.0 * M_PI * time / __orbit.orbital_period()
            );
            double secondary_x = (
                secondary_position[0] * std::cos(__arg_of_periapsis)
                -
                secondary_position[1] * std::sin(__arg_of_periapsis)
            );
            double secondary_y = (
                secondary_position[0] * std::sin(__arg_of_periapsis)
                +
                secondary_position[1] * std::cos(__arg_of_periapsis)
            );
            double center_to_secondary = secondary_position.squaredNorm();
            double position_to_secondary = (
                std::pow(position[0] - secondary_x, 2)
                +
                std::pow(position[1] - secondary_y, 2)
            );
            return (
                Core::AstroConst::G
                *
                __orbit.secondary_mass() * Core::AstroConst::solar_mass
                *
                (
                    1.0 / center_to_secondary
                    -
                    1.0 / position_to_secondary
                ) / Core::AstroConst::solar_radius
            );
        }

} //End testGravitationalPotential namespace
