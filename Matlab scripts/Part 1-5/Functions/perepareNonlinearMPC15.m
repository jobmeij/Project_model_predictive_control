function [NonLinMPC] = perepareNonlinearMPC15(Model, MPC)

yalmip('clear')

% Model data
NonLinMPC.x0 = [1 0]';
NonLinMPC.A = [1 0.1; -0.1 (1.2-0.2*NonLinMPC.x0(1)^2)];
NonLinMPC.B = [0; 0.1];               % Input is sample time Ts * u on velocity
NonLinMPC.C = [1 0];                  % Output is position
NonLinMPC.nx = 2;                     % Number of states
NonLinMPC.nu = 1;                     % Number of inputs

% MPC settings
NonLinMPC.Q = MPC.Q;
NonLinMPC.R = MPC.R;
NonLinMPC.N = MPC.N;
[~,NonLinMPC.P,~] = dlqr(NonLinMPC.A,NonLinMPC.B,NonLinMPC.Q,NonLinMPC.R);  % Compute K and P using discrete LQR
% NonLinMPC.K_lqr = -NonLinMPC.K;         
% 
% % Prediction matrices for MPC
% [NonLinMPC.Phi, NonLinMPC.Gamma] = ABN2PhiGamma(NonLinMPC.A,NonLinMPC.B,NonLinMPC.N);
% [NonLinMPC.Psi, NonLinMPC.Omega] = QRPN2PsiOmega(NonLinMPC.Q,NonLinMPC.R,NonLinMPC.P,NonLinMPC.N);
% NonLinMPC.G = 2*(NonLinMPC.Psi+NonLinMPC.Gamma'*NonLinMPC.Omega*NonLinMPC.Gamma);
% NonLinMPC.G = (NonLinMPC.G+NonLinMPC.G')/2;                                   % Making G symetric
% NonLinMPC.F = 2*NonLinMPC.Gamma'*NonLinMPC.Omega*NonLinMPC.Phi;
% % Determine invariant set
% [NonLinMPC.W,NonLinMPC.L,NonLinMPC.c] = getWLc(NonLinMPC.A,NonLinMPC.B,Model.xMax,Model.xMin,Model.uMax,Model.uMin,NonLinMPC.Gamma,NonLinMPC.Phi);

% Constraints (state and input)
NonLinMPC.umin = Model.uMin;
NonLinMPC.umax = Model.uMax;
NonLinMPC.x1min = Model.pMin;
NonLinMPC.x1max = Model.pMax;
NonLinMPC.x2min = Model.pDotMin;
NonLinMPC.x2max = Model.pDotMax;

% Creating state and input sdpvars
NonLinMPC.u = sdpvar(repmat(NonLinMPC.nu,1,NonLinMPC.N),repmat(1,1,NonLinMPC.N));
NonLinMPC.x = sdpvar(repmat(NonLinMPC.nx,1,NonLinMPC.N+1),repmat(1,1,NonLinMPC.N+1));

% Compute nonlinear MPC cost and constraints for prediction horizon
NonLinMPC.constraints = [];
NonLinMPC.cost = 0;
for k = 1:NonLinMPC.N
    NonLinMPC.cost = NonLinMPC.cost + NonLinMPC.x{k}'*NonLinMPC.Q*NonLinMPC.x{k} + NonLinMPC.u{k}'*NonLinMPC.R*NonLinMPC.u{k};
    NonLinMPC.constraints = [NonLinMPC.constraints, NonLinMPC.x{k+1} == NonlinearDiscreteModel(NonLinMPC.x{k},NonLinMPC.u{k})];
    NonLinMPC.constraints = [NonLinMPC.constraints, NonLinMPC.umin <= NonLinMPC.u{k}<= NonLinMPC.umax, NonLinMPC.x1min <=NonLinMPC.x{k+1}(1)<= NonLinMPC.x1max, NonLinMPC.x2min <=NonLinMPC.x{k+1}(2)<= NonLinMPC.x2max];
end
% Compute terminal cost and constraints
NonLinMPC.cost = NonLinMPC.cost + NonLinMPC.x{NonLinMPC.N}'*NonLinMPC.P*NonLinMPC.x{NonLinMPC.N};
NonLinMPC.constraints = [NonLinMPC.constraints, NonLinMPC.x1min <= NonLinMPC.x{k+1}(1) <= NonLinMPC.x1max, NonLinMPC.x2min <= NonLinMPC.x{k+1}(2) <= NonLinMPC.x2max];

% Inputs and outputs
NonLinMPC.parameters_in = NonLinMPC.x{1};
NonLinMPC.solutions_out = [NonLinMPC.u{:}];

% Compute control output using optimizer (nonlinear optimzation solver)
NonLinMPC.controller = optimizer(NonLinMPC.constraints,NonLinMPC.cost,[],NonLinMPC.parameters_in,NonLinMPC.solutions_out);
end

% Nonlinear discrete-time model for nonlinear MPC
function [xplus] = NonlinearDiscreteModel(x,u)
    xplus = [x(1) + 0.1*x(2); -0.1*x(1) + (1.2-0.2*x(1)^2)*x(2) + 0.1*u];
end