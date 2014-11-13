break BinarySystem.cpp:1407
commands 1
bt
print "check for lock"
print orbital_freq_mult
print spin_freq_mult
print body_index
print zone_index
print original_angmom
c
end

break BinarySystem.cpp:1419
commands 2
print "angular momentum"
print zone_ind
print zone.angular_momentum()
c
end

break BinarySystem.cpp:1429
commands 3
print "above lock fraction"
print above_lock_fraction
c
end

break BinarySystem.cpp:1436
commands 4
print spin_angmom.size()
print inclinations.size()
print periapses.size()
print spin_angmom[0]
print spin_angmom[1]
print spin_angmom[2]
print spin_angmom[3]
c
end

break RandomDiskPlanetSystem.cpp:84
commands 5
print "num locked; num_to_lock; to_lock_zone_ind; orb_mult; spin_mult"
print __num_locked_zones
print num_to_lock
print *unlocked_i
print orb_freq_mult
print spin_freq_mult
c
end

break RandomDiskPlanetSystem.cpp:103
commands 6
print "Increasing dissipation"
print __num_locked_zones
print __system->number_locked_zones()
print num_to_lock
c
end

break BinarySystem.cpp:398
commands 7
print deriv
print "matrix(0,0)="
print matrix.operator()(0,0)
print "tidal_torque_z_above[0]="
print tidal_torque_z_above[0]
print "tidal_torque_z_below[0]="
print tidal_torque_z_below[0]
c
end

break ConstPhaseLagDissipatingZone.h:103
commands 8
print "phase_lag("
print orbital_frequency_multiplier
print spin_frequency_multiplier
print forcing_frequency
print above_lock_value
c
end


run
bt
q
