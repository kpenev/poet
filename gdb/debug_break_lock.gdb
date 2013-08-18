break OrbitSolver.cpp:134
commands 1
printf "Break lock(t=%25.16g)=%25.16g (max_semi_evol=%25.16g, derivatives[0]=%25.16g)\n", age, (max_semi_evol - da_dt)/da_dt, max_semi_evol, da_dt
c
end

run
bt
q
