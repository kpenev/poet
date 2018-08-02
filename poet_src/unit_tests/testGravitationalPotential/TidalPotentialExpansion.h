/**\file
 *
 * \brief Declare an interface for evaluating the expansion of the tidal
 * potential.
 *
 * \ingroup UnitTests_group
 */

#include "../../Evolve/TidalPotentialTerms.h"

namespace testGravitationalPotential {
    ///Evaluate the tidal potential expansion.
    class TidalPotentialExpansion {
    private:
        ///\brief The coefficients of the expansion of the tidal potential.
        TidalPotentialTerms __tidal_terms;

        double
            ///The mass of the tidally perturbed object in solar masses.
            __primary_mass,

            ///The mass of the perturber object in solar masses.
            __secondary_mass,

            ///The semimajor axis of the orbit in solar radii.
            __semimajor,
            
            ///The eccentricity of the orbit
            __eccentricity,

            ///The angle from \f$ \hat{y} = \hat{S} \times \hat{L} \f$ to the
            ///direction of periapsis in radians.
            __arg_of_periapsis;
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
            
            ///See __arg_of_periapsis attribute.
            double arg_of_periapsis=Core::NaN
        ) :
            __primary_mass(primary_mass),
            __secondary_mass(secondary_mass),
            __semimajor(semimajor),
            __eccentricity(eccentricity),
            __arg_of_periapsis(arg_of_periapsis)
        {}
    }; //End TidalPotentialExpansion class.
} //End testGravitationalPotential namespace.
