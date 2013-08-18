break Common.cpp: 81

break OrbitSolver.cpp:507
commands 2
printf "Condition type=%d\n", stop_cond_type
c
end

break OrbitSolver.cpp:178
commands 3
printf "Wind saturation condition: wconv=%25.16g, wsat=%25.16g, result=%25.16g\n", wconv, wsat, (wconv-wsat)/wsat
c
end

run
bt
q
