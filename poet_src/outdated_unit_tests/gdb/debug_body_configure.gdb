break DissipatingBody.cpp:409
commands 1
print zone_index
print __orbital_frequency
print eccentricity
print orbital_angmom
print zone_angmom
print zone_inclination
print zone_periapsis
c
end

break test_DissipatingBody.cpp:80
commands 2
print body->zone(0).inclination()
print body->zone(0).angular_momentum()
print body->zone(0).moment_of_inertia(0)
c
end

break test_DissipatingBody.cpp:46
commands 3
print spin_freq_env
print inertia_env
print orbit_freq
print angmom[0]
print angmom[1]
print inclination[0]
print inclination[1]
c
end

run
q
