Implementation
==============

In order to facilitated modifying or extending the code we provide detailed
implementation documentation.

The code can be split in four parts:

 * [Stellar system](\ref StellarSystem_group):
   describing the star and the planet, including the dependence of
   non-constant quantities on the system age and other parameters.

 * [Orbit solver](\ref OrbitSolver_group): implements the steps necessary to
   calculates the orbital evolution and provides an interface for specifying
   the problem.

 * [Utilities](\ref Utilities_group): general utilities useful when defining
   the stellar system or the orbit solver.

 * [Unit tests](\ref UnitTests_group): various tests defined to make sure the
   other parts do what they claim to do.

\defgroup StellarSystem_group Stellar System
\brief Star-planet system for which the orbital evolution will be calculated.

Describes the star and the planet, including:

 * [Star](\ref Star)

   * Stellar properties (e.g. radius, moments of inertia, ...) evolving with
     age.

   * Frequency dependent tidal dissipation efficiency

   * Coupled rotational evolution of a radiative core and a convective
     envelope (for low mass stars only).

   * Loss of angular momentum due to wind

   * Having initial surface rotation locked to a proto-planetary disk
 
 * [Planet](\ref Planet)

   * Various orbital properties
   
   * Rate of tidal decay and its derivatives with respect to age and system
     or orbital parameters, as well as the torque exerted on the star due to
     the tidal dissipation.

\defgroup OrbitSolver_group Orbit Solver
\brief Calculates the orbital evolution under tides.

Solves the differential equations governing the evolution of the orbit and
the stellar spin.

\defgroup Utilities_group Utilities
\brief General utilities useful when defining the stellar system or
the orbit solver.

\defgroup UnitTests_group Unit Tests
\brief A collection of tests to verify the rest of the code.

