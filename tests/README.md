Common sense tests of POET's functionality through the Python package.
======================================================================

Stellar evolution interpolation:
--------------------------------

    [ ] Verify that higher order derivatives integrate to lower order ones.
    [ ] Verify that a collection of off-grid directly calculated tracks are
        reproduced sufficiently precisely.

Orbital evolution:
------------------
    [ ] Verify that angular momentum is conserved if wind losses are
        disabled.
    [ ] Verify that some analytically solvable evolutions are reproduced:
        [ ] Constant phase lag, non-evolving single-zone slow-spinning star.
        [ ] No dissipation no angular momentum loss two-zone star converging
            to solid body rotation
        [ ] Non-evolving star spin-down.
        [ ] Non-decaying misaligned orbit:
            [ ] [0; 90)
                [ ] S/L > 1
                [ ] S/L < 1
            [ ] [90; 180)
                [ ] S/L > 1
                [ ] S/L < 1
