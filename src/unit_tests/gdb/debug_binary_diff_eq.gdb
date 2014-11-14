break DissipatingZone.cpp:400 if deriv==0
commands 1
print "mp, m"
print mp
print m
print "worb, wspin, wtide"
print orbital_frequency
print __spin_frequency
print tidal_frequency
print "phase lags"
print mod_phase_lag_above
print mod_phase_lag_below
print "locks"
print __lock
print __other_lock
c
end

break test_BinarySystem.cpp:45
commands 2
print mp
print m
print system.orbital_frequency()
print system.quantity(static_cast<Quantity>(FIRST_ZONE_ANGVEL+zone_ind))
print forcing_freq
print system.__lags[zone_ind].operator()(m, mp)
c
end

run
bt
q
