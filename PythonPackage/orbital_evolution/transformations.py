"""For now just define transformation from log10(Q) to phase lag."""

from numpy import pi, log10

#lgQ is the most intuitive name for log10(Q).
#pylint: disable=invalid-name
#Names reflect functions are each others' inverses
#pylint: disable=redefined-outer-name
def phase_lag(lgQ):
    """Return the phase lag corresponding to the given Q value."""

    return 15.0 / (16.0 * pi * 10.0**lgQ)

def lgQ(phase_lag):
    """Return the log10(Q) value corresponding to the given phase lag."""

    return log10(15.0 / (16.0 * pi * phase_lag))
#pylint: enable=invalid-name
#pylint: enable=redefined-outer-name
