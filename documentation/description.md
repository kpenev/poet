Description {#mainpage}
===========

<strong>P</strong>lanetary <strong>O</strong>rbital <strong>E</strong>volution
due to <strong>T</strong>ides (POET) is a tool for simulating the evolution of
single stars, binary stars, and star--planet systems under the combined effects
of:

 * Tides in one or both of the objects

 * Age dependent stellar structure

 * Stars losing angular momentum to wind

 * The internal redistribution of angular momentum between the surface and the interior of stars

POET is capable of evolution calculations for orbits ranging from circular to
highly eccentric, and for fully flexible prescription for the tidal dissipation
efficiency.

Briefly, stars in POET are split into a number of zones each of which
experiences:

  * Boundaries evolving both in mass and radius

  * Its own tidal coupling to the companion

  * Exchange of angular momentum with its neighboring zones, both due to
    "friction" and due to the exchange of material as the zone boundaries
    evolve.

In addition, the surface zone loses angular momentum due to winds, modeled
after, [Stauffer and Hartmann
1987](https://ui.adsabs.harvard.edu/abs/1987ApJ...318..337S), [Kawaler
1988](https://ui.adsabs.harvard.edu/abs/1988ApJ...333..236K), and [Barnes and
Sofia 1996](https://ui.adsabs.harvard.edu/abs/1996ApJ...462..746B), as

\f[
    \dot{L} = K_w
              \omega \min(\omega, \omega_{sat})^2
              \sqrt{(R_\star/R_\odot)(M_\odot/M_\star)}
\f]

where \f$K_w\f$ parameterizes the strength of the wind, \f$\omega_{sat}\f$ is a
saturation frequency required by observations, and \f$R_\star\f$ and
\f$M_\star\f$ are the mass and radius (evolving) of the star. The coupling
torque between neighboring zones is such that it would cause exponential decay
of the differential rotation on a timescale \f$\tau_c\f$ (separately defined for
each pair of zones) if everything else was turned off.

Stellar evolution is incorporated using interpolation within a grid of stellar
properties as a function of mass, metallicity and age. 

At the moment, POET comes with a pre--computed grid of stellar evolutions for
interpolation for \f$0.4 M_\odot \leq M_\star < 1.4 M_\odot\f$, which was
generated using the publicly available stellar evolution code MESA ([Paxton et.
al. 2011](http://adsabs.harvard.edu/abs/2011ApJS..192....3P)).

Tidal evolution uses an extended version of the formalism descrsibed in [Lai
2012](http://adsabs.harvard.edu/abs/2012MNRAS.423..486L), which was extended to
allow for eccentric orbits, where a series expansion in eccentricity is used,
up to a dynamically adjusted maximum order.

For more details see:

  * [Installation Instructions](\ref Installation)

  * [Python interface](\ref PythonInterface)

  * [C++ interface](\ref CPPInterface)

  * [A detailed derivation of the equations specifying the evolution](\ref InclinationEccentricity).

  * [Acknowledgements](\ref Acknowledgements)

  * If you use POET in your publications, and/or presentations, please cite
    [Penev, Zhang and Jackson
    2014](http://adsabs.harvard.edu/abs/2014arXiv1405.1050P)
    and the
    [Astrophysics Source Code Library entry](https://ui.adsabs.harvard.edu/abs/2014ascl.soft08005P/abstract)
