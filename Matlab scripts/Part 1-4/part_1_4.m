% MPC graded homework assignment
% Part 1.4 -  Quasi LPV MPC 
% Job Meijer            - 1268155
% Marcel van Wensveen   - 1253085


% Init script
clear all
close all
close_system    % Close all open simulink models in memory
clc

% Add functions directory
addpath('./functions/')

%% Init variables for simulation quasiLPV
syms p;
p =sym('p');
Ad = [1, 0.1;-0.1, 1.2-p];
QuasiLPV.Ad = matlabFunction(subs(Ad,p,p), 'vars', p);
QuasiLPV.Bd = [0; 0.1];
QuasiLPV.Cd = eye(2);
QuasiLPV.Dd = [0];
QuasiLPV.p_eq = @(x) 0.2*x.^2;

QuasiLPV.Q = [0.1 0; 0 1];
QuasiLPV.R = 0.01;
QuasiLPV.N = 10;

% Simulation parameters
QuasiLPV.x0 =  [1; 0];
QuasiLPV.SimT = 20;
QuasiLPV.Tau = 0.1;

% Constraints
QuasiLPV.x_low = [-2; -5];
QuasiLPV.x_high = [2; 5];
QuasiLPV.u_low = [-1];
QuasiLPV.u_high = [1];

%% Set params for lin MPC and nonlin MPC

% Set simulation settings
SimSettings.x0 = QuasiLPV.x0;                   % Initial condition
SimSettings.Ts = QuasiLPV.Tau;                  % Sample time [s]
SimSettings.T0 = 0;                             % Start time [s]
SimSettings.Tend = QuasiLPV.SimT;               % End time [s]
SimSettings.t = SimSettings.T0:SimSettings.Ts:SimSettings.Tend;  % Time vector
SimSettings.DTplot = false;                     % Plotting trajectories for DT model simulation
SimSettings.CTplot = false;                     % Plotting trajectories for CT model simulation
SimSettings.LinMPCplot = false;
SimSettings.NonLinMPCplot = false;

% Model settings
Model.uMin = QuasiLPV.u_low;                    % Input min constraint
Model.uMax = QuasiLPV.u_high;                   % Input max constraint
Model.pMin = QuasiLPV.x_low(1);                 % State p min constraint
Model.pMax = QuasiLPV.x_high(1);                % State p max constraint
Model.pDotMin = QuasiLPV.x_low(2);              % State p' min constraint
Model.pDotMax = QuasiLPV.x_high(2);             % State p' max constraint
Model.xMin = [Model.pDotMin; Model.pMin];       % Combining
Model.xMax = [Model.pDotMax; Model.pMax];       % Combining

% Linearized CT ss model:       
Model.Act = [0 1;-1 0];                         % Removed P^2(t) from equation (lin. @ [1; 0])
Model.Bct = [0; 1];                             % Input at state 2 (pddot)
Model.Cct = [1 0];                              % Only output p(t) (state 1)
Model.Dct = 0;

Model.LinCtSys = ss(Model.Act,Model.Bct,Model.Cct,Model.Dct);
Model.LinDTSys = c2d(Model.LinCtSys,SimSettings.Ts);

Model.Alin = Model.LinDTSys.A;        
Model.Blin = Model.LinDTSys.B;            
Model.Clin = Model.LinDTSys.C;             
Model.Dlin = Model.LinDTSys.D;

% MPC controller settings (for closed loop linear and nonlinear MPC simulations)
MPC.N = QuasiLPV.N;                                         % Prediction horizon
MPC.Q = QuasiLPV.Q;                                         % prediction cost which makes linear MPC feasable (1.2)        
MPC.R = QuasiLPV.R;                                         % Input cost
[MPC.K,MPC.P,~] = dlqr(Model.Alin,Model.Blin,MPC.Q,MPC.R);  % Compute K and P using discrete LQR
MPC.K_lqr = -MPC.K;                                         % Invert K
% Prediction matrices for MPC
[MPC.Phi, MPC.Gamma] = ABN2PhiGamma(Model.Alin,Model.Blin,MPC.N);
[MPC.Psi, MPC.Omega] = QRPN2PsiOmega(MPC.Q,MPC.R,MPC.P,MPC.N);
MPC.G = 2*(MPC.Psi+MPC.Gamma'*MPC.Omega*MPC.Gamma);
MPC.G = (MPC.G+MPC.G')/2;                                   % Making G symetric
MPC.F = 2*MPC.Gamma'*MPC.Omega*MPC.Phi;
% Determine invariant set
[MPC.W,MPC.L,MPC.c] = getWLc(Model.Alin,Model.Blin,Model.xMax,Model.xMin,Model.uMax,Model.uMin,MPC.Gamma,MPC.Phi);


%% Calculate invariant, control adm., terminal set and cost matrices
[QuasiLPV] = Calculate_sets_QuasiLPV(QuasiLPV);

%% Simulate with quasiLPV
QuasiLPV = Simulate_QuasiLPV(QuasiLPV, false);

%% Simulation of linear MPC in closed loop with Simulink model
[LinMPC] = runLinearMPC(SimSettings,Model,MPC);

%% Simulation of nonlinear MPC in closed loop with Simulink model
[NonLinMPC] = runNonlinearMPC(SimSettings,Model,MPC);

%% Comparing quasiLPV, linear MPC and nonlinear MPC controllers in closed loop
compareMPCcontrollers_14(LinMPC,NonLinMPC, QuasiLPV);

