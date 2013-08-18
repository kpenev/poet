break OrbitSolver.cpp:1141
commands 1
set $tstart=t
set $tmax=max_next_t
set $trewind_max=NaN
set $trewound_actual=NaN
c
end

break OrbitSolver.cpp:1144
commands 2
set $tend=t
printf "GSL step from %25.16g to %25.16g (limit %25.16g).\n", $tstart, $tend, $tmax
c
end

break OrbitSolver.cpp:1148
commands 3
set $orig_stop=stop
c
end

break OrbitSolver.cpp:1150
commands 4
set $trewind_max=stop.stop_age()
set $trewound_actual=t
c
end

break OrbitSolver.cpp:1159
commands 5
printf "Rejected=%d step from %25.16g to %25.16g, limited to %25.16g. Original stop=\n", step_rejected, $tstart, $tend, $tmax
print $orig_stop
printf "Updated stop=\n"
print stop
printf "Rewound from %25.16g to %25.16g (max requested %25.16g)\n", $tend, $trewound_actual, $trewind_max
c
end

break OrbitSolver.cpp:1165
commands 6
print "GSL error"
bt
end

break OrbitSolver.cpp:1168
commands 7
print "step accepted"
c
end

break OrbitSolver.cpp:1182
commands 8
printf "Evolution stopped at t=%25.16g for mode change.", t
c
end

run
bt
q
