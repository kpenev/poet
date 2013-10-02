break OrbitSolver.cpp:1042
commands 1
	printf "t=%25.16e, max_next_t=%25.16e\n", t, max_next_t
	c
end

set args -M 0.9000000000000001 -m 10.00000000000001 -p 4.611111111111112 --p-disk 1.4 -K 0.155 --core-env-coupling-timescale=12 --t-disk 2.5 --low-mass-wind-sat-w=2.454 -r 0.714 --precision 5 --serialized-stellar-evol interp_state_data_phs4 --tmax 10

run
bt
q
