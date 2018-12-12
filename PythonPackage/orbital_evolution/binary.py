"""Define a class for binaries which can be evolved."""

from ctypes import c_int, c_bool, c_double
import numpy

from orbital_evolution.star_interface import EvolvingStar
from orbital_evolution.planet_interface import LockedPlanet
from orbital_evolution.evolve_interface import library
from basic_utils import\
    Structure,\
    calc_semimajor,\
    calc_orbital_frequency,\
    calc_orbital_angular_momentum

def get_evolution_quantities(secondary_is_star):
    """
    Return the list of quantities for a binary system.

    Args:
        secondary_is_star: True iff the evolution was for a binary star (vs a
        planet-star) system.

    Returns:
        A list of the evolution quantities tracked for an evolution of the
        given system type.
    """

    evolution_quantities = ['age', 'semimajor', 'eccentricity']

    star_float_quantities = ['envelope_inclination',
                             'core_inclination',
                             'envelope_periapsis',
                             'core_periapsis',
                             'envelope_angmom',
                             'core_angmom']

    for quantity in star_float_quantities:
        if secondary_is_star:
            evolution_quantities.append('primary_' + quantity)
        else:
            evolution_quantities.append(quantity)


    if secondary_is_star:
        for quantity in star_float_quantities:
            evolution_quantities.append('secondary_' + quantity)

    evolution_quantities.append('evolution_mode')

    if secondary_is_star:
        evolution_quantities.extend(['primary_wind_saturation',
                                     'secondary_wind_saturation'])
    else:
        evolution_quantities.append('wind_saturation')

    return evolution_quantities

#Two of the "attributes are actually methods
#pylint: disable=too-many-instance-attributes
class Binary:
    """A class for binaries POET can evolve."""

    evolution_modes = ['LOCKED_SURFACE_SPIN',
                       'BINARY',
                       'SINGLE',
                       'TABULATION']
    _evolution_mode_ids = {
        mode: c_int.in_dll(library, mode + '_EVOL_MODE').value
        for mode in evolution_modes
    }

    @staticmethod
    def evolution_quantity_c_type(quantity):
        """Return the ctypes type of the given evolution quantity."""

        if quantity == 'evolution_mode':
            return c_int
        if quantity.endswith('wind_saturation'):
            return c_bool
        return c_double

    def __init__(self,
                 primary,
                 secondary,
                 *,
                 disk_lock_frequency,
                 disk_dissipation_age,
                 initial_semimajor=None,
                 initial_orbital_period=None,
                 initial_eccentricity=0.0,
                 initial_inclination=0.0,
                 secondary_formation_age=None):
        """
        Create a binary system out of two bodies.

        Args:
            - primary:
                The first body in the system. Assumed to always be there, so
                for a star-planet system this should be the star.

            - secondary:
                The second body in the system, initially may not be there and
                later may be engulfed by the first body.

            - disk_lock_frequency:
                Frequency of the surface spin of the primary when disk is
                present in rad/day.

            - disk_dissipation_age:
                Age when disk dissipates in Gyrs.

            - initial_semimajor:
                The semimajor axis of the orbit at which the secondary forms
                in solar radii. If omitted, initial_orbital_period must be
                specified.

            - initial_orbital_period:
                Alternative to specifying the initial semimajor axis.

            - initial_eccentricity:
                The eccentricity of the orbit at which the secondary forms.

            - initial_inclination:
                Inclination between surface zone of primary and initial orbit
                in radians.

            - secondary_formation_age:
                Age when the secondary forms.

        Returns: None
        """

        assert isinstance(primary, EvolvingStar)
        self.primary = primary
        self.secondary = secondary
        if initial_semimajor is None:
            initial_semimajor = self.semimajor(initial_orbital_period)
        if secondary_formation_age is None:
            secondary_formation_age = disk_dissipation_age

        self.evolution_quantities = get_evolution_quantities(
            isinstance(secondary, EvolvingStar)
        )
        if isinstance(secondary, LockedPlanet):
            c_create_func = library.create_star_planet_system
            self._c_get_evolution_func = library.get_star_planet_evolution
            self._c_get_final_state = library.get_star_planet_final_state
        else:
            assert isinstance(secondary, EvolvingStar)
            c_create_func = library.create_star_star_system
            self._c_get_evolution_func = library.get_star_star_evolution
            self._c_get_final_state = library.get_star_star_final_state

        self.c_binary = c_create_func(
            primary.c_body,
            secondary.c_body,
            initial_semimajor,
            initial_eccentricity,
            initial_inclination,
            disk_lock_frequency,
            disk_dissipation_age,
            secondary_formation_age
        )

        self.num_evolution_steps = 0
        self.c_solver = None

    def delete(self):
        """Destroy the binary created at construction."""

        library.destroy_binary(self.c_binary)
        if hasattr(self, 'c_solver'):
            library.destroy_solver(self.c_solver)

    def configure(self,
                  *,
                  age,
                  semimajor,
                  eccentricity,
                  spin_angmom,
                  inclination,
                  periapsis,
                  evolution_mode):
        """
        Set the current state (orbit) of a system.

        Args:
            - age:
                The age to set the system to.

            - semimajor:
                The semimajor axis of the orbit in solar radii.

            - eccentricity:
                The eccentricity of the orbit.

            - spin_angmom:
                The spin angular momenta of the zones of the bodies (body 1
                first, outermost zone to innermost, followed by body 2).

            - inclination:
                The inclinations of the zones of the bodies (same order as
                spin_angmom). The surface zone inclination must be omitted
                for single body systems.

            - periapsis:
                The arguments of periapsis of the zones of the bodies (same
                order as spin_angmom, but not including the surface zone of
                the first body).

            - evolution_mode:
                The evolution mode to assume. Must be one of the constants
                defined.

        Returns: None
        """

        library.configure_system(self.c_binary,
                                 age,
                                 semimajor,
                                 eccentricity,
                                 spin_angmom,
                                 inclination,
                                 periapsis,
                                 self._evolution_mode_ids[evolution_mode])

    def evolve(self,
               final_age,
               max_time_step,
               precision,
               required_ages):
        """
        Evolve the system forward from its current state.

        Args:
            - final_age:
                The age at which to stop the evolution in Gyrs. The starting
                age must be already set for the system through configure.

            - max_time_step:
                The maximum size of the time step allowed in Gyrs.

            - precision:
                The precision to require of the solution.

            - required_ages:
                Ages at which the evolution must stop precisely.

        Returnns: None
        """

        #The point is to check if previous call defined the member
        #pylint: disable=access-member-before-definition
        if hasattr(self, 'c_solver'):
            library.destroy_solver(self.c_solver)
        #pylint: enable=access-member-before-definition
        self.c_solver = library.evolve_system(
            self.c_binary,
            final_age,
            max_time_step,
            precision,
            required_ages,
            (0 if required_ages is None else required_ages.size)
        )
        self.num_evolution_steps = library.num_evolution_steps(self.c_solver)

    def get_evolution(self, quantities=None):
        """
        Return the last calculated evolution.

        Args:
            - quantities:    An iterable of quantities to read the evolution of.
                The evolution of omitted quantities can still be obtained later
                by subsequent calls to this method. The allowed entries are in
                the star_star_evolution_quantities for binary star systems or in
                star_planet_evolution_quantities for a star-planet system. If
                None, it defaults to the full list of quantities for the given
                system.

        Returns:
            Sturture:
                A structure with mebers named the same way as the input list of
                quantities containing the values of the corresponding quantity
                at each evolution step. The order is always in increasing age.
        """

        result = Structure()

        if quantities is None:
            quantities = self.evolution_quantities

        for quantity_name in quantities:
            setattr(
                result,
                quantity_name,
                numpy.empty(
                    self.num_evolution_steps,
                    dtype=self.evolution_quantity_c_type(quantity_name)
                )
            )

        get_evol_args = [self.c_solver,
                         self.c_binary,
                         self.primary.c_body]
        if isinstance(self.secondary, EvolvingStar):
            get_evol_args.append(self.secondary.c_body)
        get_evol_args.extend([getattr(result, quantity, None)
                              for quantity in self.evolution_quantities])
        self._c_get_evolution_func(*get_evol_args)
        return result

    def final_state(self):
        """Return the final evolution state of a system (all quantities)."""

        result = Structure()
        library_final_state = [self.evolution_quantity_c_type(q)()
                               for q in self.evolution_quantities]
        if isinstance(self.secondary, LockedPlanet):
            self._c_get_final_state(self.c_solver,
                                    self.c_binary,
                                    self.primary.c_body,
                                    *library_final_state)
        else:
            self._c_get_final_state(self.c_solver,
                                    self.c_binary,
                                    self.primary.c_body,
                                    self.secondary.c_body,
                                    *library_final_state)
        for quantity, library_value in zip(self.evolution_quantities,
                                           library_final_state):
            setattr(result, quantity, library_value.value)
        return result

    def orbital_frequency(self, semimajor):
        """
        The orbital frequency of the system for the given semimajor axis.

        Args:
            - semimajor:
                The semimajor axis at which the system's orbital period is
                required in solar radii.

        Returns:
            The orbital period in days if the two bodies of this system are
            in an orbit with the given semimajor axis.
        """

        return calc_orbital_frequency(self.primary.mass,
                                      self.secondary.mass,
                                      semimajor)

    def orbital_period(self, semimajor):
        """
        The orbital period of the system for the given semimajor axis.

        Args:
            - semimajor:
                The semimajor axis at which the system's orbital period is
                required in solar radii.

        Returns:
            The orbital period in days if the two bodies of this system are
            in an orbit with the given semimajor axis.
        """

        return 2.0 * numpy.pi / self.orbital_frequency(semimajor)

    def semimajor(self, orbital_period):
        """
        The semimajor axis of the system for the given orbital period.

        Args:
            - orbital_period:
                The orbital period at which the system's orbital period is
                required in days.

        Returns:
            The semimajor axis in solar radii if the two bodies of this
            system are in an orbit with the given period.
        """

        return calc_semimajor(self.primary.mass,
                              self.secondary.mass,
                              orbital_period)

    def orbital_angular_momentum(self, semimajor, eccentricity):
        """
        The orbital agular momentum for the given semimajor/eccentricity.

        Args:
            - semimajor:
                The semimajor axis of the system.

            - eccentricity:
                The orbital eccentricity.

        Returns:
            The orbital angular momentum if the two bodies are in an orbit
            with the given semimajor axis and eccentricity in solar units.
        """

        return calc_orbital_angular_momentum(self.primary.mass,
                                             self.secondary.mass,
                                             semimajor,
                                             eccentricity)
#pylint: enable=too-many-instance-attributes
