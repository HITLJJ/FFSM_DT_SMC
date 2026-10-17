function B = input_matrix()
%INPUT_MATRIX Map the two internal joint torques into five generalized forces.
B = [zeros(3,2); eye(2)];
end
