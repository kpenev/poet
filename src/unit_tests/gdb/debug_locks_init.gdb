break BinarySystem.cpp:389 if deriv==0
commands 1
print "lock, nontidal_torque, tidal_torque_above, tidal_torque_below, tidal power diff:"
print nontidal_torque[0]
print tidal_torque_z_above[0]
print tidal_torque_z_below[0]
print tidal_power_difference[0]
print "rhs="
print rhs[0]
c
end

break BinarySystem.cpp:572
commands 2
print "global orbit torque="
print global_orbit_torque[0]
print global_orbit_torque[1]
print global_orbit_torque[2]
print "zone_orbit_torque="
print zone_orbit_torque[0]
print zone_orbit_torque[1]
print zone_orbit_torque[2]
print "total_zone_torque="
print total_zone_torque[0]
print total_zone_torque[1]
print total_zone_torque[2]
c
end

break BinarySystem.cpp:1429
commands 3
print "below lock multiplier("
print orbital_freq_mult
print spin_freq_mult
print body_index
print zone_index
print ")="
print 2.0*above_lock_fraction-1.0
c
end

break test_BinarySystem.cpp:707
commands 4
print "Expected below lock coef("
print locked_zone_ind
print ")="
print lock_coef
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

break test_BinarySystem.cpp:184
commands 5
print "Expected orbit frequency rate="
print -1.5*expected_diff_eq[0]*worb/a
print "Orbutal angular momentum="
print orbital_angmom
c
end

break DissipatingBody.cpp:350 if deriv==0
commands 6
print "orbit torque below:"
print __orbit_torque[0][0]
print __orbit_torque[0][1]
print __orbit_torque[0][2]
print "orbit torque correction"
print correction[0]
print correction[1]
print correction[2]
c
end

run
bt
q
