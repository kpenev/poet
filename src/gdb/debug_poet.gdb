#break YRECEnvelope::configure
#commands 1
#print "Configuring YREC envelope"
#c
#end
#
#break YRECCore::configure
#commands 2
#print "Configuring YREC Core"
#c
#end
#
#break LockedPlanet::configure
#commands 3
#print "Configuring locked planet"
#c
#end
#
#break DissipatingZone.cpp:359
#commands 4
#print __spin_frequency
#print __orbital_frequency
#print __lock
#print __other_lock
#if isnan(inclination)
#	bt
#end
#c
#end
#
#break DissipatingZone.cpp:566
#commands 5
#print __lock
#c
#end

#break WindSaturationCondition.cpp:25
#commands 1
#print evol_mode
#print __primary
#print __body.number_zones()
#print __other_body.number_zones()
#print angmom_index
#print surf_angmom_deriv
#print __body.zone(0).moment_of_inertia(1)
#print wsurf
#print __body.zone(0).moment_of_inertia(0)
#print __saturation_freq
#c
#end

break DissipatingBody.cpp:458 if deriv==0 && isnan(result[2])
commands 1
print result
print zone_index
c
end

break DissipatingBody.cpp:459 if deriv==0 && isnan(result[1])
commands 2
print result
print zone_index
c
end

break DissipatingBody.cpp:462 if deriv==0 && isnan(result[1])
commands 3
print result
print zone_index
c
end

break DissipatingBody.cpp:467 if deriv==0 && isnan(result[1])
commands 4
print result
print zone_index
c
end

run
bt
q
