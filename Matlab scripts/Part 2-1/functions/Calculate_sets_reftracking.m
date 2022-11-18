function [refMPC] = Calculate_sets_reftracking(refMPC)

disp('Start calculating sets and linear cost matrices for linear MPC')

% Calculate terminal set and controller
[refMPC.P, ~, ~] = dare(refMPC.Ad, refMPC.Bd, refMPC.Q, refMPC.R);

% Create MPC matrices
disp('Create cost matrices')
[refMPC] = create_predic_matric_reftracking(refMPC);
[refMPC] = create_cost_matric_reftracking(refMPC);
[refMPC] = create_gain_matric_reftracking(refMPC);
[refMPC] = Calculate_reftracking_matrices(refMPC);
[refMPC] = Construct_MPC_constraints_reftracking(refMPC);

disp('Finished calculating sets and linear cost matrices for linear MPC')
