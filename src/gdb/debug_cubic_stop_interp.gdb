break OrbitSolver.cpp:481
commands 1
print "Extremum without derivatives"
c
end

break OrbitSolver.cpp:490
commands 2
print "Extremum with derivatives"
c
end

break OrbitSolver.cpp:446
commands 3
print "Quadratic extremum from 3 points"
c
end

break OrbitSolver.cpp:463
commands 4
print "Cubic extremum from 4 points"
c
end

break OrbitSolver.cpp:637
commands 5
print t
print stop
c
end

break Common.cpp:204
commands 6
printf "First cubic extremum at x=%25.15g\n", extremum_x
c
end

break Common.cpp:206
commands 7
printf "Final cubic extremum at x=%25.16g\n", extremum_x
c
end

break Common.cpp:206
condition 8 (extremum_x<require_range_low || extremum_x>require_range_high)
commands 8
printf "extremum x=%25.16g\n", extremum_x
end

run
bt
q
