break DissipatingZone.cpp:291
commands 1
print __spin_frequency/__orbital_frequency
print __lock
print __other_lock
c
end

break test_DissipatingBody.cpp:23
commands 2
print a
print m_env
print m_core
print r_env
print r_core
print m_other
print inertia_env
print inertia_core
print spin_freq_env
print spin_freq_core
print orbit_freq
c
end

run
bt
q
