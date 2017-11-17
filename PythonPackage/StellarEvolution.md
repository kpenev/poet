The Python Stellar Evolution Module {#StellarEvolutionReadMe}
===================================

A module for interpolating among existing stellar evolution tracks. Since
generating the interpolation is quite slow, POET supports "serializing"
(storing as a file on disk) a computed interpolation. The
`stellar_evolution.manager` sub-module keeps track of exsiting
interpolations, reusing them whenever appropriate, so this is the suggested
interface.

### Basic usage :

The user requests an interpolator from the manager by specifying either the
interpolator name (using the `StellarEvolutionManager.get_interpolator_by_name`
method) or by specifying a configuration for the interpolator. In the former
case, it is an error to request an intorpolator with a name not previously
defined (POET comes with a pre-defined interpolator named `default`). In the
latter case, if the user specifies a name to assign to the new interpolator, the
manager will automatically generate an interpolator if one does not exist.

    from poet.stellar_evolution.manager import StellarEvolutionManager
    import numpy

    %Create a manager
    manager = StellarEvolutionManager(<path to serialized interpolators>)

    %Find and use the default interpolator
    interpolator = manager.get_interpolator_by_name('default')

    % The convective zone moment of inertia for a 0.73 solar mass star with
    % [Fe/H] = -0.13
    iconv = interpolator('ICONV', 0.73, -0.13)

    % Evaluate the moment of inertia at a list of ages.
    iconv_values = iconv(numpy.linspace(0.05, 10.0, 100))

### Custom interpolation using the package tracks:

If for some reason the standard interpolation is not optimal for your
purposes, it is possible to modify the parameters for how each quantity's age
interpolation is performed (the subsequent interpolation in mass and
metallicity is always done using a bi-cubic spline):

    % construct an interpolator using only a subset of the tracks, as well
    % as custom smoothing and number of nodes.
    interpolator = manager.get_interpolator(
        nodes = dict(RADIUS = <int>,    %number of nodes for the radius age
                                        %interpolation
                     ICONV  = <int>,    %convective moment of inertia nodes
                     LUM    = <int>,    %luminosity nodes
                     IRAD   = <int>,    %radiative moment of inertia nodes
                     MRAD   = <int>,    %radiative mass nodes
                     RRAD   = <int>),   %tachocline radius nodes
        smoothing = dict(RADIUS = <float>,  %smoothing parameter for the
                                            %radius age interpolation
                         ICONV  = <float>,  %convective moment of inertia
                                            %smoothig
                         LUM    = <float>,  %luminosity smoothing
                         IRAD   = <float>,  %radiative moment of inertia
                                            %smoothing
                         MRAD   = <float>,  %radiative mass smoothing
                         RRAD   = <float>), %tachocline radius smoothing
        vs_log_age = dict(RADIUS = <bool>,  %should the independent argument
                                            %for radius interpolation be
                                            %log(age) instead of age?
                         ICONV  = <bool>,   %Same as RADIUS but for the
                                            %convective moment of inertia.
                         LUM    = <bool>,   %same as RADIUS but for luminosity.
                         IRAD   = <float>,  %same as RADIUS but for radiative
                                            %moment of inertia.
                         MRAD   = <float>,  %same as RADIUS but for mass in the
                                            %radiative core.
                         RRAD   = <float>), %same as RADIUS but for radius in
        log_quantity = dict(RADIUS = <bool>,  %should the log(radius) be used
                                              %in the interpolation instead of
                                              %radius?
                            ICONV  = <bool>,  %Same as RADIUS but for the
                                              %convective moment of inertia.
                            LUM    = <bool>,  %same as RADIUS but for
                                              %luminosity.
                            IRAD   = <float>, %same as RADIUS but for radiative
                                              %moment of inertia.
                            MRAD   = <float>, %same as RADIUS but for mass in
                                              %the radiative core.
                            RRAD   = <float>),%same as RADIUS but for radius in
                                              %the radiative core.
        track_fnames = [...], %List of filenames containing stellar evolution
                              %tracks to interpolate among.
        masses = [<float>, ...],        %list of stellar masses to use
        feh = [<float>, ...]            %list of stellar [Fe/H] values to use
        model_suite = <str>,            %The name of a preiously registered
                                        %suite of tracks to use in the
                                        %interpolation. POET comes with a suite
                                        %named 'MESA'.
        new_interp_name = <str> %A human-readable name to assign should an
                                %interpolator matching the other arguments
                                %does not exist.
        num_threads = <int> %The number of threads to use when deriving the
                            %interpolation.
    )

As explained above, all subsequent calls to `get_interpolator` using the same
parameters, will re-use the interpolator generated by the first call. Note
the `new_interp_name` argument. If that argument is not supplied, and an
interpolator matching the specified configuration is not found, the result is
None.

### Solving for stellar mass and age which match other known properties ###

Another useful task performed by the [Stellar Evolution](#StellarEvolutionReadMe)
module is to solver for the mass and age of a star given [Fe/H] and any two of:
\f$T_{eff}\f$, \f$\log_{10}(g)\f$, \f$L_\star\f$, \f$\rho_\star\f$. This is done
using the interpolators change_varables method. For example:

    from poet.stellar_evolution.manager import StellarEvolutionManager
    import numpy

    %Create a manager
    manager = StellarEvolutionManager(<path to serialized interpolators>)

    %Find and use the default interpolator
    interpolator = manager.get_interpolator_by_name('default')

    solutions = manager.change_variables(
        feh = 0.0,
        teff = 6000,
        lum = 1.0
    )
    for mass, age in solutions :
        ....

### Custom stellar evolution tracks

In is also possible to construct interpolations based on a set of tracks not
packaged with POET. The tracks must be formatted as comma separated values,
with the first line containing the names of the columns.  All quantities
required by POET must be tabulated. See the set of tracks shipped with POET
for an example. 

The first step is to register the tracks with the interpolation manager. This
can be done track by track:
    
    %Register a single stellar evolution track.
    manager.register_track(<filename>,  %Name of the file contaning the track
                           mass,        %The stellar mass of the track (Msun)
                           metallicity, %[Fe/H] of the track.
                           model_suite) %Name of a track collection to assign
                                        %this track to.

If the names of the files can be parsed to extract the mass and metallicity,
there is a shortcut:

    %Register an entire collection of tracks at once with stellar mass and
    %metallicity encoded in the filename.
    manager.register_track_collection(
        [<filename>,...],     %List of track files
        <regular expression>, %A compiled python re which defines groups
                              %names 'MASS' and either 'Z' or '[Fe/H]'
                              %The package tracks match:
                              %'M(?P<MASS>[0-9.E+-]+)_Z(?P<Z>[0-9.E+-]+).csv'
        model_suite           %Name of a track collection to assign the
                              %tracks to
    )

After registering the tracks, new interpolators can be generated and
registered with the manager, either like we did above using the custom
tracks, but this time we specify our new model suite.

    % construct an interpolator using only a subset of the tracks, as well
    % as custom smoothing and number of nodes.
    interpolator = manager.get_interpolator(
        masses = [<float>, ...],        %list of stellar masses to use
        metallicities = [<float>, ...], %list of stellar metallicities to use
        model_suite = <str>,            %the model suite defined above
        nodes = dict(RADIUS = <int>,    %number of nodes for the radius age
                                        %interpolation
                     ICONV  = <int>,    %convective moment of inertia nodes
                     LUM    = <int>,      %luminosity nodes
                     IRAD   = <int>,    %radiative moment of inertia nodes
                     MRAD   = <int>,    %radiative mass nodes
                     RRAD   = <int>),   %tachocline radius nodes
        smoothing = dict(RADIUS = <float>,  %smoothing parameter for the
                                            %radius age interpolation
                         ICONV  = <float>,  %convective moment of inertia
                                            %smoothig
                         LUM    = <float>,  %luminosity smoothing
                         IRAD   = <float>,  %radiative moment of inertia
                                            %smoothing
                         MRAD   = <float>,  %radiative mass smoothing
                         RRAD   = <float>), %tachocline radius smoothing
        new_interp_name = <str> %A human-readable name to assign should an
                                %interpolator matching the other arguments
                                %not exist.
    )

This only works if the new suite has only one track per mass-metallicity
combination. If for example the new tracks again use `model_suite = MESA`,
then which tracks to use must be specified by name:

    % construct an interpolator using a custom set of tracks, as well
    % as custom smoothing and number of nodes.
    interpolator = manager.get_interpolator(
        track_fnames = [<path>, ...],   %list of filenames containing
                                        %the stellar evolution tracks
        nodes = <same as above>         
        smoothing = <same as above>
        new_interp_name = <same as above>
    )
