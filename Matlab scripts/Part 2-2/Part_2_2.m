% MPC graded homework assignment
% Part 2.2 -  MPC design reference tracking with disturbance rejection
% Job Meijer            - 1268155
% Marcel van Wensveen   - 1253085

% Goal: Design a linear MPC controller with disturbance model to track the
% given reference with zero offset. Plotting the estimated disturbance from
% the observer and the actual disturbance, as well as the closed-loop
% output trajectories over the desired reference.

% Init
clear all; close all; clc;

% Add functions directory
addpath('./functions/')

% System
refMPC.Ad = [0.9719 0.0155; 0.2097 0.9705];
refMPC.Bd = [0.2097; 0.3263];
refMPC.BdDist = [0.5; 0.2];
refMPC.Cd = [1 0];
refMPC.Dd = [0];
refMPC.g1 = [0 0];

% Create disturbance array
refMPC.dk = [0.1.*ones(600,1); 0.5.*ones(600,1); -0.5.*ones(601,1)];

% Cost
refMPC.Q = [100 0; 0 0.001];
refMPC.R = 1;
refMPC.N = 18; % was 18

% Simulation parameters
refMPC.x0 =  [8; 0];
refMPC.SimT = 1800;
refMPC.Tau = 1;

% Constraints
refMPC.x_low = [-15; -100];
refMPC.x_high = [30; 100];
refMPC.u_low = [-25];
refMPC.u_high = [25];

% Disturbance estimation and rejection
refMPC.x_dist_est = zeros(2,length(0:refMPC.Tau:refMPC.SimT));
refMPC.d_dist_est = zeros(1,length(0:refMPC.Tau:refMPC.SimT));
refMPC.Q_dist_est = eye(2);
refMPC.R_dist_est = 1;
refMPC.Lx_dist = [0.5; 0.5];
refMPC.Ld_dist = 0.1;
refMPC.x_ss = zeros(2,length(0:refMPC.Tau:refMPC.SimT));
refMPC.u_ss = zeros(1,length(0:refMPC.Tau:refMPC.SimT));

% Reference
refMPC.ref = [10*ones(600,1);-10*ones(600,1);10*ones((600+refMPC.N+1),1)];


%% Calculate invariant, control adm., terminal set and cost matrices
refMPC = Calculate_sets_reftracking(refMPC);


%% Simulate with linear discrete time model
refMPC = Simulate_refMPC_dist(refMPC, true);

