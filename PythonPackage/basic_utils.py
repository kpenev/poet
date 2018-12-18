"""A collection of orbital functions and a class for storing data."""

from math import pi
from astropy import units, constants

#Class intended as data representation.
#pylint: disable=too-few-public-methods
class Structure:
    """An empty class used only to hold user defined attributes."""

    def __init__(self, **initial_attributes):
        """Create a class with (optionally) initial attributes."""

        for attribute_name, attribute_value in initial_attributes.items():
            setattr(self, attribute_name, attribute_value)

    def format(self, prefix=''):
        """Generate a tree-like representation of self."""

        result = ''
        for attr_name in dir(self):
            if attr_name[0] != '_' and attr_name != 'format':
                attribute = getattr(self, attr_name)
                if isinstance(attribute, Structure):
                    result += (prefix
                               +
                               '|-'
                               +
                               attr_name
                               +
                               '\n'
                               +
                               attribute.format(prefix + '| '))
                else:
                    result += (prefix
                               +
                               '|-'
                               +
                               attr_name
                               +
                               ': '
                               +
                               str(attribute)
                               +
                               '\n')
        return result
#pylint: enable=too-few-public-methods

#In this case m1 and m2 make perfect sense as argument names.
#pylint: disable=invalid-name
def calc_semimajor(m1, m2, orbital_period):
    """
    Return the semimajor axis for a binary with the given orbit.

    Args:
        - m1:
            Mass of the first body in the system

        - m2:
            Mass of the second body in the system.

        - orbital_period:
            The orbital period of the binary in days.

    Returns:
        The semimajor axis in solar radii.
    """

    #Pylint false positive, does not detect units and constants members.
    #pylint: disable=no-member
    return (
        (
            constants.G * (m1 + m2) * units.M_sun
            *
            (orbital_period * units.day)**2
            /
            (4.0 * pi**2)
        )**(1.0 / 3.0)
    ).to(units.R_sun).value
    #pylint: enable=no-member
#pylint: enable=invalid-name

#In this case m1 and m2 make perfect sense as argument names.
#pylint: disable=invalid-name
def calc_orbital_frequency(m1, m2, semimajor):
    """
    Return the orbital for a binary with the given orbit.

    Args:
        - m1:
            Mass of the first body in the system

        - m2:
            Mass of the second body in the system.

        - semimajor:
            The semimajor axis of the binary in solar radii.

    Returns:
        The orbital angular velocity in rad / day.
    """

    #Pylint false positive, does not detect units and constants members.
    #pylint: disable=no-member
    return (
        (
            constants.G * (m1 + m2) * units.M_sun
            /
            (semimajor * units.R_sun)**3
        )**0.5
    ).to(1 / units.day).value
    #pylint: enable=no-member
#pylint: enable=invalid-name

#In this case m1 and m2 make perfect sense as argument names.
#pylint: disable=invalid-name
def calc_orbital_angular_momentum(m1, m2, semimajor, eccentricity):
    """
    Return the angular momentum of the given orbit.

    Args:
        - m1:
            Mass of the first body in the system

        - m2:
            Mass of the second body in the system.

        - semimajor:
            The semimajor axis of the binary in solar radii.

        - eccentricity:
            The eccentricity of the orbit.

    Returns:
        The orbital angular momentum in solar units for a binary with the
        given bodies in the given orbit.
    """

    #Pylint false positive, does not detect units and constants members.
    #pylint: disable=no-member
    return (
        m1 * m2 / (m1 + m2) * semimajor**2
        *
        calc_orbital_frequency(m1, m2, semimajor)
        *
        (1.0 - eccentricity**2)**0.5
    )
    #pylint: enable=no-member
#pylint: enable=invalid-name
