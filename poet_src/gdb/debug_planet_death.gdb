break OrbitSolver.cpp:576
commands 1
printf "age=%25.16g, extremum_x=%25.16g, stop_cond_value=%25.16g, extremum.y()=%25.16g, stop_cond_history.back()[cond_ind]=%25.16g\n", age, extremum.__x, stop_cond_value, extremum.__y, stop_cond_history.back()[cond_ind]
c
end

break OrbitSolver.cpp:645
commands 2
print stop
c
end

break OrbitSolver.cpp:532
commands 3
printf "Zero crossing from p1=(%25.16g, %25.16g), (%25.16g, %25.16g) and d1=%25.16g, d2=%25.16g\n", prev_age, prev_stop_cond, age, stop_condition_value, prev_stop_deriv, stop_condition_derivative
c
end

break OrbitSolver.cpp:482
commands 4
print result
c
end

break Common.cpp:173
commands 5
printf "Cubic extremum: extremum_x=%25.16g, extremum_y=%25.16g, a=%25.16g, b=%25.16g, c=%25.16g, d=%25.16g\n", x, *extremum_y, a, b, c, d
c
end

break Common.cpp:177
commands 6
printf "Linear extremum: extremum_x=%25.16g, extremum_y=%25.16g\n", extremum, *extremum_y
c
end

run
q
