"""For now just define transformation from log10(Q) to phase lag."""

from math import pi

#This is the most intuitive name for log10(Q).
#pylint: disable=invalid-name
def phase_lag(lgQ):
    """Return the phase lag corresponding to the given Q value."""

    return 15.0 / (16.0 * pi * 10.0**lgQ)
#pylint: enable=invalid-name
