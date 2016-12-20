A set of modules exposing POET's functionality in python.
=========================================================
 - [Stellar Evolution](#stellar_evolution)

<h2 id = 'stellar_evolution'>Stellar Evolution (stellar_evolution)</h2>
-------------------------------------

A module for interpolating among existing stellar evolution tracks. Since
generating the interpolation is quite slow, POET supports "serializing"
(storing as a file on disk) a computed interpolation. The
`stellar_evolution.manager` module keeps track of exsiting interpolations,
reusing them whenever appropriate, so this is the suggested interface. 

### Usage example :

    from poet.stellar_evolution.manager import StellarEvolutionManager
    import numpy

    %Create a manager
    manager = StellarEvolutionManager()

    % construct an interpolator using the default configuration.
    interpolator = manager.get_interpolator() 

    % The convective zone moment of inertia for a 0.73 solar mass star with 
    % [Fe/H] = -0.13
    iconv = interpolator('ICONV', 0.73, -0.13)

    % Evaluate the moment of inertia at a list of ages.
    iconv_values = iconv(numpy.linspace(0.05, 10.0, 100))

### Customizing the interpolation:
If for some reason the standard interpolation is not optimal for your
purposes, it is possible to modify the collection of tracks to use, as well
as the parameters for how each quantity's age interpolation is performed:

    % construct an interpolator using only a subset of the tracks, as well
    % as custom smoothing and number of nodes.
    interpolator = manager.get_interpolator(
        masses = [<float>, ...],        %list of stellar masses to use
        metallicities = [<float>, ...]  %list of stellar metallicities to use
        nodes = dict(RADIUS = <int>,    %number of nodes for the radius age
                                        %interpolation
                     ICONV  = <int>,    %convective moment of inertia nodes
                     LUM    = <int>,      %luminosity nodes
                     IRAD   = <int>,    %radiative moment of inertia nodes
                     MRAD   = <int>,    %radiative mass nodes
                     RRAD   = <int>),   %tachocline radius nodes
        smoothing = dict(RADIUS = <int>,    %smoothing parameter for the
                                            %radius age interpolation
                         ICONV  = <int>,    %convective moment of inertia
                                            %smoothig
                         LUM    = <int>,    %luminosity smoothing
                         IRAD   = <int>,    %radiative moment of inertia
                                            %smoothing
                         MRAD   = <int>,    %radiative mass smoothing
                         RRAD   = <int>),   %tachocline radius smoothing
    )

In is also possible to construct interpolations based on a set of tracks not
packaged with POET. At this point the tracks must be formatted as comma
separated values, with the first line containing the names of the columns.
All quantities required by POET must be tabulated. See the set of tracks
shipped with POET for an example.

    % construct an interpolator using a custom set of tracks, as well
    % as custom smoothing and number of nodes.
    interpolator = manager.get_interpolator(
        track_fnames = [<path>, ...],   %list of filenames containing
                                        %the stellar evolution tracks
        track_suite = <str>,            %Name to assign to this set of tracks
                                        %e.g. YY or YREC. Note that the
                                        %tracks shipped with POET have
                                        %track_suite = 'MESA'
        nodes = <same as above>         
        smoothing = <same as above>
    )
