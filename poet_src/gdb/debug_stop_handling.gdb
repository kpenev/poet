break OrbitSolver.cpp: 477
commands 1
printf "--> gsl_odeiv2_evolve_apply(t=%+23.16e, old_t=%+23.16e, step=%+23.16e, orbit(%d==%d, 0th=%+23.16e), max_next_t=%+23.16e, evolution_mode=%d\n", t, old_t, step_size, orbit.size(), nargs, orbit._M_data[0], max_next_t, evolution_mode
c
end

break OrbitSolver.cpp: 261
commands 2
printf "--> diff_eq(age=%+23.16e, mode=%d, wind=%d)=%+23.16e, retcode=%d\n", age, evol_mode, wind_state, orbital_derivatives[0], result
c
end

break OrbitSolver.cpp: 297
commands 3
printf "--> Jac(age=%+23.16e, orbit=%+23.16e, mode=%d, wind=%d)=%+23.16e, %+23.16e, retcode=%d\n", age, orbital_parameters[0], evol_mode, wind_state, param_derivs[0], age_derivs[0], result
end

break OrbitSolver.cpp: 491
commands 4
printf "--> Stop information: stop_reason=%d, crossing=%d, stop_age=%+23.16e, stop=%+23.16e, t=%+23.16e, old_t=%+23.16e, max_next_t=%+23.16e", stop_info.stop_reason, stop_info.is_crossing, stop_info.stop_age, stop, t, old_t, max_next_t
c
end

break /Volumes/SVNDiskImage/tidal_orbital_evolution/trunk/Common.cpp:86
commands 5
printf "--> x0=%+23.16e, x1=%+23.16e, x=%+23.16e, y0=%+23.16e, y1=%+23.16e, dy0=%+23.16e, dy1=%+23.16e", x0, x1, x, y0, y1, dy0, dy1
c
end

break /Volumes/SVNDiskImage/tidal_orbital_evolution/trunk/Common.cpp:93
commands 6
printf "--> x0=%+23.16e, x1=%+23.16e, extremum=%+23.16e, y0=%+23.16e, y1=%+23.16e, dy0=%+23.16e, dy1=%+23.16e", x0, x1, x, y0, y1, dy0, dy1
c
end

break OrbitSolver.cpp: 493
commands 7
printf "--> going back to t=%+23.16e", old_t
c
end

break OrbitSolver.cpp: 509
commands 8
printf "--> accepting t=%+23.16e", t
c
end

break OrbitSolver.cpp: 595
commands 9
printf "--> next_evol_mode(t=%+23.16e, old_mode=%d, condition=%d, condition_value=%+23.16e, stopped_before=%d\n", age, evolution_mode, condition_type, condition_value, stopped_before
c
end

break OrbitSolver.cpp: 907
commands 10
printf "--> Changing mode from %d to %d\n", old_evolution_mode, evolution_mode
c
end

delete 2
delete 3

run
q
