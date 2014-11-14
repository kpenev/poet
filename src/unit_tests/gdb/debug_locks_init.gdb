break BinarySystem.cpp:1424
commands 1
bt
print "check for lock"
print orbital_freq_mult
print spin_freq_mult
print body_index
print zone_index
print original_angmom
print body.zone(zone_index).__lock
print __semimajor
print __body1.mass()
print __body2.mass()
c
end

break BinarySystem.cpp:1429
commands 2
print "above lock fraction"
print above_lock_fraction
c
end

break RandomDiskPlanetSystem.cpp:100
commands 3
print "Increasing dissipation"
print __num_locked_zones
print __system->number_locked_zones()
print num_to_lock
c
end

#break DissipatingZone.cpp:340 if this->__lock.__lock_direction==0
#commands 4
#print "Configuring locked zone"
#print orbital_frequency
#print eccentricity
#print orbital_angmom
#print spin_angmom
#print inclination
#print periapsis
#print __orbital_angmom
#print __orbital_frequency
#print __spin_frequency
#print __angular_momentum
#c
#end

#break DissipatingZone.cpp:416 if deriv_ind==0 && mp==2 && m==1 && __lock.__lock_direction==0
#commands 5
#print __lock
#print mp
#print m
#print term_torque_z
#print mod_phase_lag_below
#print mod_phase_lag_above
#print __torque_z[deriv_ind]
#print __torque_z[deriv_ind+1]
#c
#end

run
bt
q
