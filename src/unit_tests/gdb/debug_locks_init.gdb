break DissipatingZone.cpp:173
commands 1
print __spin_frequency/__orbital_frequency
print __lock
print __other_lock
print __other_lock.lock_direction()*__other_lock.spin_frequency_multiplier()*spin_frequency()
print __other_lock.lock_direction()*__other_lock.orbital_frequency_multiplier()*__orbital_frequency
c
end

run
bt
q
