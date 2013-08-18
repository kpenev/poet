break OrbitSolver.cpp:750
commands 1
printf "pre1_cond=%25.16g, pre2_cond=%25.16g, post_cond=%25.16g", pre1_cond, pre2_cond, post_cond
c
end

break OrbitSolver.cpp:618
commands 2
printf "go_back=%d, *first_stop_cond[%d]=%25.16g\n", go_back, cond_ind, (*first_stop_cond)[cond_ind]
c
end

run
bt
q
