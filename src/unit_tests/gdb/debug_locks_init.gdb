break BinarySystem.cpp:1429
commands 1
print "below lock multiplier("
print orbital_freq_mult
print spin_freq_mult
print body_index
print zone_index
print ")="
print 2.0*above_lock_fraction-1.0
c
end

break test_BinarySystem.cpp:741
commands 2
print "Expected below lock coef("
print locked_zone_ind
print ")="
print lock_coef
print "expected rhs, expected matrix, expected orbital and spin terms"
print numer-denom/2.0
print denom
print denom1
print denom2
print "locked_spin_dir="
print locked_spin_dir[0]
print locked_spin_dir[1]
print locked_spin_dir.norm()
print "locked spin angmom="
print locked_inertia*system_maker.orbital_frequency()
print "locked, nonlocked torques, orbit torques, orbit_coef for locked:"
print locked_tidal_torque[0]
print locked_tidal_torque[1]
print zone_torques[locked_zone_ind][0]
print zone_torques[locked_zone_ind][1]
print orbit_torque[0]
print orbit_torque[1]
print orbit_coef
print "locked spin dir="
print locked_spin_dir[0]
print locked_spin_dir[1]
print "Expected spin frequency rate="
print (zone_torques[locked_zone_ind][0]*locked_spin_dir[0]+zone_torques[locked_zone_ind][1]*locked_spin_dir[1])/locked_inertia
c
end

break test_BinarySystem.cpp:186
commands 3
print "Expected orbit frequency rate="
print -1.5*expected_diff_eq[0]*worb/a
print "Orbutal angular momentum="
print orbital_angmom
c
end

break BinarySystem.cpp:411 if deriv==0
commands 4
print "Computing above lock fraction from rhs, matrix:"
print rhs[0]
print matrix.operator()(0,0)
c
end

break BinarySystem.cpp:393 if deriv==0
commands 5
print "Orbit term="
print 1.5*tidal_power_difference[i]/__orbital_energy
c
end

break BinarySystem.cpp:398 if deriv==0
commands 6
print "Spin term="
print (tidal_torque_z_above[i]-tidal_torque_z_below[i])/zone.angular_momentum()
print "Ztorque above, below, angmom"
print tidal_torque_z_above[i]
print tidal_torque_z_below[i]
print zone.angular_momentum()
c
end

run
bt
q
