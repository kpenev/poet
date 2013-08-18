break OrbitSolver.cpp: 452
condition 1 stop>0
commands 1
print evolution_mode
print wind_state
print stop_reason
printf "stop_cond(%+23.16e)=%+23.16e, old_stop_cond(%+23.16e)=%+23.16e, max_next_t=%+23.16e, stop=%+23.16e\n", t, *(&stop_condition_value), old_t, old_stop_condition_value, max_next_t, stop
end

break OrbitSolver.cpp: 417
condition 2 stop_cond.type(i)==SYNCHRONIZED
commands 2
print evolution_mode
printf "sync_cond(%+23.16e)=%+23.16e\n", t, *(&current_stop_cond[i])
c
end

break OrbitSolver.cpp: 114
commands 3
printf "age=%+23.16e, dwconv_dt=%+23.16e, max_semi_evol=%+23.16e, derivatives[0]=%+23.16e, stop_cond=%+23.16e\n", age, dwconv_dt, max_semi_evol, *(&derivatives[0]), (max_semi_evol-derivatives[0])/derivatives[0]
end

break OrbitSolver.cpp: 386
commands 4
print stop_cond.type(i)
printf "Zeroing old_stop_cond[%d](%+23.16e)=%+23.16e\n", i, t, old_stop_cond[i]
c
end

break OrbitSolver.cpp: 390
commands 4
print stop_cond.type(i)
printf "Changing sign of old_stop_cond[%d](%+23.16e)=%+23.16e\n", i, t, old_stop_cond[i]
c
end
