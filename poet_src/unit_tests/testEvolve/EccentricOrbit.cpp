/**\file
 *
 * \brief Define some of the methods of EccentricyOrbit.
 *
 * \ingroup UnitTests_group
 */

#include "EccentricOrbit.h"

double eccentric_anomaly_equation(double eccentric_anomaly,
                                  void *orbital_phase_eccentricity)
{
    double orbital_phase = reinterpret_cast<double*>(
        orbital_phase_eccentricity
    )[0];
    double eccentricity = reinterpret_cast<double*>(
        orbital_phase_eccentricity
    )[1];
    return (eccentric_anomaly
            -
            eccentricity * std::sin(eccentric_anomaly)
            -
            orbital_phase);
}

namespace Evolve {
    double EccentricOrbit::reduced_mass() const
    {
        return (
            (__primary_mass * __secondary_mass)
            /
            (__primary_mass + __secondary_mass)
        );
    }

    double EccentricOrbit::eccentric_anomaly(double orbital_phase) const
    {
        const gsl_root_fsolver_type *FsolverType;
        gsl_root_fsolver *solver;
        double root = 0;
        double lower_bound = 0.0, upper_bound = 2.0 * M_PI;
        gsl_function Function;
        orbital_phase -= 2.0 * M_PI * std::floor(orbital_phase / (2.0 * M_PI));
        double params[] = {orbital_phase, __eccentricity};

        Function.function = &eccentric_anomaly_equation;
        Function.params = params;

        FsolverType = gsl_root_fsolver_brent;
        solver = gsl_root_fsolver_alloc (FsolverType);
        gsl_root_fsolver_set(solver, &Function, lower_bound, upper_bound);

#ifdef VERBOSE_DEBUG
        std::cerr << std::setw(25) << "Iteartion"
                  << std::setw(25) << "LowerBound"
                  << std::setw(25) << "root"
                  << std::setw(25) << "UpperBound"
                  << std::endl;
#endif

        int status = GSL_CONTINUE;
        for(unsigned iter = 0; status == GSL_CONTINUE; ++iter) {
            status = gsl_root_fsolver_iterate(solver);
            root = gsl_root_fsolver_root(solver);
            lower_bound = gsl_root_fsolver_x_lower(solver);
            upper_bound = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(lower_bound,
                                            upper_bound,
                                            1e-12,
                                            1e-12);
#ifdef VERBOSE_DEBUG
            std::cerr << std::setw(25) << iter
                      << std::setw(25) << lower_bound
                      << std::setw(25) << root
                      << std::setw(25) << upper_bound
                      << std::endl;
#endif
        }

        gsl_root_fsolver_free(solver);

        assert(status == GSL_SUCCESS);

        return root;
    }

    Eigen::Vector3d EccentricOrbit::secondary_position(
        double orbital_phase
    ) const
    {
        double current_eccentric_anomaly = eccentric_anomaly(orbital_phase);
        return Eigen::Vector3d(
            __semimajor * (std::cos(current_eccentric_anomaly)
                           -
                           __eccentricity),
            __semimajor * (sqrt(1.0 - std::pow(__eccentricity, 2))
                           *
                           std::sin(current_eccentric_anomaly)),
            0.0
        );
    }

    double EccentricOrbit::orbital_angmom() const
    {
        return (
            reduced_mass() * Core::AstroConst::solar_mass
            *
            std::pow(__semimajor * Core::AstroConst::solar_radius, 2)
            *
            (2.0 * M_PI / (orbital_period() * Core::AstroConst::day))
            *
            std::sqrt(1.0 - std::pow(__eccentricity, 2))
        );
    }

    double EccentricOrbit::orbital_energy() const
    {
        return (
            -Core::AstroConst::G
            *
            __primary_mass * Core::AstroConst::solar_mass
            *
            __secondary_mass * Core::AstroConst::solar_mass
            /
            (2.0 * __semimajor * Core::AstroConst::solar_radius)
        );
    }

    double EccentricOrbit::orbital_period() const
    {
        return std::sqrt(
            4.0 * M_PI * M_PI
            *
            std::pow(__semimajor * Core::AstroConst::solar_radius, 3)
            /
            (
                Core::AstroConst::G
                *
                (__primary_mass + __secondary_mass)
                *
                Core::AstroConst::solar_mass
            )
        ) / Core::AstroConst::day;
    }
}
