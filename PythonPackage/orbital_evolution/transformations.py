from math import pi
from astropy import constants, units

def phase_lag(lgQ) :
    """Return the phase lag corresponding to the given Q value."""

    return 15.0 / (16.0 * pi * 10.0**lgQ);
