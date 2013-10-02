break OrbitSolver.cpp: 475
commands 1
printf "DISCARDING t=%+23.16e -> %+23.16e, mode=%3d, stop_reason=%3d", old_t, t, evolution_mode, stop_reason
c
end

break OrbitSolver.cpp: 483
commands 2
printf "ACCEPTING  t=%+23.16e -> %+23.16e, mode=%3d, stop_reason=%3d", old_t, t, evolution_mode, stop_reason
c
end

run
bt
q
