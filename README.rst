****************************************
Planetary Orbital Evolution due to Tides
****************************************

Python wrapped C++ code for calculating the orbital evolution of two bodies
under mutually raised tides.

`Documentation on GitHubPages <https://kpenev.github.io/poet/>`_

Three-zone tidal dissipation
============================

This branch lets POET set tidal dissipation separately in the convective
envelope (zone 0), radiative shell (zone 1), and convective core (zone 2).
Rebuild the native POET libraries before using the Python interface. Use
``convert_history_to_track`` to make paired ``radiative_shell`` and
``convective_core`` tracks from the same MESA history and age range. Pass
their interpolators to ``EvolvingStar`` as ``interpolator`` and
``convective_core_interpolator``, then configure each zone with
``set_dissipation``.

See ``scripts/run_equal_mass_three_zone.py`` for an example that takes a
MESA history and an eccentricity coefficient database.

Potentially Interesting Articles
================================

  * ` A self-consistent weak friction model for the tidal evolution of circumbinary planets  <https://ui.adsabs.harvard.edu/abs/2019arXiv190605195Z/abstract>`_
  
Tidal Love Number
=================

  * `Kopal (1959) "Close Binary Systems" <https://archive.org/details/closebinarysyste00kopa/page/27/mode/1up>`_ 

  * `Claret (2019) "Updating the theoretical tidal evolution constants: Apsidal motion and the moment of inertia" <https://ui.adsabs.harvard.edu/abs/2019A%26A...628A..29C/abstract>`_
  
  * `Cook (1980) "Interiors of the Planets" <https://github.com/kpenev/poet/blob/master/documentation/Interiors%20of%20the%20Planets%20(Cambridge%20Planetary%20Science%20Old)%20by%20A.%20H.%20Cook%20(z-lib.org).pdf>`_
