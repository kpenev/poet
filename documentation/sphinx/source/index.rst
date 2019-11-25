.. (P)lanetary (O)rbital (E)volution due to (T)ides documentation master file, created by
   sphinx-quickstart on Fri Nov 22 09:57:25 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


*******************************************************
(P)lanetary (O)rbital (E)volution due to (T)ides (POET)
*******************************************************

**P**\ lanetary **O**\ rbital **E**\ volution due to **T**\ ides (POET) is a
tool for simulating the evolution of single stars, binary stars, and
star--planet systems under the combined effects of:

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
after, `Stauffer and Hartmann 1987
<https://ui.adsabs.harvard.edu/abs/1987ApJ...318..337S>`_, `Kawaler 1988
<https://ui.adsabs.harvard.edu/abs/1988ApJ...333..236K>`_, and `Barnes and Sofia
1996 <https://ui.adsabs.harvard.edu/abs/1996ApJ...462..746B>`_, as

.. math::

    \dot{L} = K_w
              \omega \min(\omega, \omega_{sat})^2
              \sqrt{(R_\star/R_\odot)(M_\odot/M_\star)}

where :math:`K_w` parameterizes the strength of the wind, :math:`\omega_{sat}`
is a saturation frequency required by observations, and :math:`R_\star` and
:math:`M_\star` are the mass and radius (evolving) of the star. The coupling
torque between neighboring zones is such that it would cause exponential decay
of the differential rotation on a timescale :math:`\tau_c` (separately defined
for each pair of zones) if everything else was turned off.

Stellar evolution is incorporated using interpolation within a grid of stellar
properties as a function of mass, metallicity and age. 

At the moment, POET comes with a pre--computed grid of stellar evolutions for
interpolation for :math:`0.4 M_\odot \leq M_\star < 1.4 M_\odot`\ , which was
generated using the publicly available stellar evolution code MESA (`Paxton et.
al. 2011 <http://adsabs.harvard.edu/abs/2011ApJS..192....3P>`_).

Tidal evolution uses an extended version of the formalism descrsibed in `Lai
2012 <http://adsabs.harvard.edu/abs/2012MNRAS.423..486L>`_, which was extended
to allow for eccentric orbits, where a series expansion in eccentricity is used,
up to a dynamically adjusted maximum order.

For more details see:

.. toctree::
    :maxdepth: 1

    installation

    Python Interface <_implementation/modules>

    The C++ interface <cpp_api/library_root>

    A detailed derivation of the equations specifying the evolution <inclination_eccentricity_equations>

    acknowledgements

Citing POET
===========

If you use POET in your publications, and/or presentations, please cite `Penev,
Zhang and Jackson 2014 <http://adsabs.harvard.edu/abs/2014arXiv1405.1050P>`_ and
the `Astrophysics Source Code Library entry
<https://ui.adsabs.harvard.edu/abs/2014ascl.soft08005P/abstract>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
