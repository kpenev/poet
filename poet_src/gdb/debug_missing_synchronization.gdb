break OrbitSolver.cpp:1088
commands 1
printf "Mode change requested at t=%25.16g, reason=", last_age
print stop_reason
printf "stop condition value=%25.16g, stopped_before=", stop_condition_value
print stopped_before
printf "Old evolution mode="
print evolution_mode
printf "Mode change is forbidden: %d", no_evol_mode_change
c
end

break OrbitSolver.cpp:1094
commands 2
printf "New evolution mode="
print evolution_mode
c
end

run
q
