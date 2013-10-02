break poet.cpp:606
commands 1
	print real_parameters[InCol::MAX_STEP]
	c
end

break poet.cpp:550
commands 2
	print real_parameters[InCol::MAX_STEP]
	c
end

break poet.cpp:716
commands 3
	print real_parameters[InCol::MAX_STEP]
	c
end

break poet.cpp:686
commands 4
	print real_parameters[InCol::MAX_STEP]
	c
end

break poet.cpp:340
commands 5
	print __direct_value_options[InCol::MAX_STEP]->dval[0]
	print __defaults[InCol::MAX_STEP]
	c
end

set args -M 0.9000000000000001 -m 10.00000000000001 -p 4.611111111111112 --p-disk 1.4 -K 0.155 --core-env-coupling-timescale=12 --t-disk 2.5 --low-mass-wind-sat-w=2.454 -r 0.714 --precision 5 --serialized-stellar-evol interp_state_data_phs4 --tmax 10

run
bt
q
