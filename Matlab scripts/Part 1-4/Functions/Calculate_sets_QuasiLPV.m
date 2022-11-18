function [QuasiLPV] = Calculate_sets_QuasiLPV(QuasiLPV)

disp('Start calculating sets and linear cost matrices for Quasi LPV')

% Calculate terminal set and controller
[QuasiLPV.P, ~, K] = dare(QuasiLPV.Ad(QuasiLPV.p_eq(QuasiLPV.x0(1))), QuasiLPV.Bd, QuasiLPV.Q, QuasiLPV.R);
QuasiLPV.K = -K;
QuasiLPV.Acl = QuasiLPV.Ad(QuasiLPV.p_eq(QuasiLPV.x0(1)))+QuasiLPV.Bd*QuasiLPV.K;
% 
% Calculate the invariant set
disp('Calculate invariant set')
QuasiLPV.LTImodel = LTISystem('A', QuasiLPV.Acl); %Acl
QuasiLPV.Xset = Polyhedron([-eye(size(QuasiLPV.x_low,1)); eye(size(QuasiLPV.x_low,1))], [-QuasiLPV.x_low;QuasiLPV.x_high]);
QuasiLPV.XUset = Polyhedron([-eye(size(QuasiLPV.u_low,1));eye(size(QuasiLPV.u_low,1))]*QuasiLPV.K,[-QuasiLPV.u_low; QuasiLPV.u_high]) & QuasiLPV.Xset;
QuasiLPV.InvSetXU = QuasiLPV.LTImodel.invariantSet('X',QuasiLPV.XUset);

% Create MPC matrices
disp('Create cost matrices')
QuasiLPV.Mn = QuasiLPV.InvSetXU.A;
QuasiLPV.bn = QuasiLPV.InvSetXU.b;
[QuasiLPV] = create_predic_matric_quasi(QuasiLPV);
[QuasiLPV] = create_cost_matric_quasi(QuasiLPV);
[QuasiLPV] = create_gain_matric_quasi(QuasiLPV);
[QuasiLPV] = Construct_MPC_constraints_with_terminal_quasi(QuasiLPV);

disp('Finished calculating sets and linear cost matrices for linear MPC')
