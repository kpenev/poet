****************************************
Planetary Orbital Evolution due to Tides
****************************************

Python wrapped C++ code for calculating the orbital evolution of two bodies
under mutually raised tides.

`Documentation on GitHubPages <https://kpenev.github.io/poet/>`_

POET for massive stars
======================

This branch updates a stellar zone boundary algorithm for massive stars.
At each MESA model step the algorithm identifies the central convective
core, the intervening radiative shell, and the outer convective envelope
from the mesh's mixing types. It starts at the center and finds the first
unmixed or anonymously mixed cell outside the convective core. From the
surface it scans blocks of up to 100 mesh cells for a block with no
convective mixing, then locates the first unmixed or anonymously mixed
cell below the last convective cell above that block. The MESA history
records both boundary locations and the moments of inertia of the three
regions.

When converting that history to a POET track, the default ``normal`` mode
chooses a two-zone boundary at each time step. If the convective core
has a larger moment of inertia than the convective envelope, the
radiative shell is grouped with the envelope. Otherwise, the radiative
shell is grouped with the core. Use ``rc_env`` or ``rc_core`` to keep the
shell assigned to the envelope or core, respectively, for the whole
track. The two physical boundaries come from the MESA
``run_star_extras.f`` template. The grouping is implemented by
``convert_history_to_track`` in
``PythonPackage/stellar_evolution/MIST_realtime/build_model.py``.

Start with the linked POET build documentation. To generate new MIST2
tracks, supply a MIST2 work directory grid through ``mist2_base_dir`` and
configure the MESA runtime for your system. Bring your own generated
tracks and eccentricity coefficient database for orbital runs that use
those inputs.

Potentially Interesting Articles
================================

  * ` A self-consistent weak friction model for the tidal evolution of circumbinary planets  <https://ui.adsabs.harvard.edu/abs/2019arXiv190605195Z/abstract>`_
  
Tidal Love Number
=================

  * `Kopal (1959) "Close Binary Systems" <https://archive.org/details/closebinarysyste00kopa/page/27/mode/1up>`_ 

  * `Claret (2019) "Updating the theoretical tidal evolution constants: Apsidal motion and the moment of inertia" <https://ui.adsabs.harvard.edu/abs/2019A%26A...628A..29C/abstract>`_
  
  * `Cook (1980) "Interiors of the Planets" <https://github.com/kpenev/poet/blob/master/documentation/Interiors%20of%20the%20Planets%20(Cambridge%20Planetary%20Science%20Old)%20by%20A.%20H.%20Cook%20(z-lib.org).pdf>`_
