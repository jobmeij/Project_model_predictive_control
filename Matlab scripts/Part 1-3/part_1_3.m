% MPC graded homework assignment
% Part 1.3 -  Nonlinear MPC design
% Job Meijer            - 1268155
% Marcel van Wensveen   - 1253085

close all;
clear all;  
close_system    % Close all open simulink models in memory
clc;


% Add functions directory
addpath('./functions/')

% Set simulation settings
SimSettings.x0 = [1 0]';            % Initial condition
SimSettings.Ts = 0.1;               % Sample time [s]
SimSettings.T0 = 0;                 % Start time [s]
SimSettings.Tend = 20;              % End time [s]
SimSettings.t = SimSettings.T0:SimSettings.Ts:SimSettings.Tend;  % Time vector
SimSettings.DTplot = false;         % Plotting trajectories for DT model simulation
SimSettings.CTplot = true;          % Plotting trajectories for CT model simulation
SimSettings.LinMPCplot = false;
SimSettings.NonLinMPCplot = false;

% Model settings
Model.uMin = -1;                    % Input min constraint
Model.uMax = 1;                     % Input max constraint
Model.pMin = -2;                    % State p min constraint
Model.pMax = 2;                     % State p max constraint
Model.pDotMin = -5;                 % State p' min constraint
Model.pDotMax = 5;                  % State p' max constraint
Model.xMin = [Model.pDotMin; Model.pMin];       % Combining
Model.xMax = [Model.pDotMax; Model.pMax];       % Combining

% Linearized CT ss model:
Model.Act = [0 1;-1 0];             % Removed P^2(t) from equation (lin. @ [1; 0])
Model.Bct = [0; 1];                 % Input at state 2 (pddot)
Model.Cct = [1 0];                  % Only output p(t) (state 1)
Model.Dct = 0;

Model.LinCtSys = ss(Model.Act,Model.Bct,Model.Cct,Model.Dct);
Model.LinDTSys = c2d(Model.LinCtSys,SimSettings.Ts);

Model.Alin = Model.LinDTSys.A;        
Model.Blin = Model.LinDTSys.B;            
Model.Clin = Model.LinDTSys.C;             
Model.Dlin = Model.LinDTSys.D;

% MPC controller settings (for closed loop linear and nonlinear MPC simulations)
MPC.N = 50;                                                 % Prediction horizon
MPC.Q = [0.1 0; 0 1];                                       % prediction cost which makes linear MPC feasable (1.2)        
MPC.R = 0.01;                                               % Input cost
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


%% Run discrete-time nonlinear model
disp('Simulating discrete time nonlinear model')
[DTsim] = runDiscreteTimeModel(SimSettings);     

%% Run continuous-time nonlinear (Simulink) model
disp('Simulating continious time nonlinear model')
[CTsim] = runContinuousTimeModel(SimSettings, Model, MPC);

%% Comparing discrete- and continuous time simulations
compareDtCtSimulations(DTsim,CTsim,SimSettings);

%% Simulation of linear MPC in closed loop with Simulink model
disp('Simulating linear MPC with nonlinear model')
[LinMPC] = runLinearMPC(SimSettings,Model,MPC);

%% Simulation of nonlinear MPC in closed loop with Simulink model
disp('Simulating nonlinear MPC with nonlinear model')
[NonLinMPC] = runNonlinearMPC(SimSettings,Model,MPC);

%% Comparing linear and nonlinear MPC controllers in closed loop
compareMPCcontrollers_13(SimSettings,LinMPC,NonLinMPC);
