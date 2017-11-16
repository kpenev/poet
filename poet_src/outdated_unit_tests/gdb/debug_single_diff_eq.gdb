break test_BinarySystem.cpp:168
commands 1
print expected_diff_eq[0]
print expected_diff_eq[1]
print expected_diff_eq[2]
print expected_diff_eq[3]
print angmom_loss
c
end

break DissipatingBody.cpp:460
commands 2
print angular_momentum_loss(0)
c
end

break BinarySystem.cpp:103
commands 3
print zone_index
print torque[0]
print torque[1]
print torque[2]
c
end

run
bt
q
