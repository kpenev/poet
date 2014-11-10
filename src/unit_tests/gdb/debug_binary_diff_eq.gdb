break DissipatingZone.cpp:493
commands 1
print orbit_x_torque
print orbit_z_torque
print __orbital_angmom
print (orbit_x_torque*cos_inc-orbit_z_torque*sin_inc)/__orbital_angmom
print zone_x_torque
print zone_x_torque/__angular_momentum
print result
c
end

break test_BinarySystem.cpp:364
commands 2
print x_torque
print zone_rotation
c
end

break test_BinarySystem.cpp:367
commands 3
print orbit_rotation
print orbital_angmom
print worb
c
end

run
bt
q
