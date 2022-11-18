% MPC graded homework assignment
% Part 1.5 -  calculate feasable set for quasi LPV, linMPC and non-lin MPC 
% Job Meijer            - 1268155
% Marcel van Wensveen   - 1253085


% Init script
clear all
close all
close_system    % Close all open simulink models in memory
clc

% Add functions directory
addpath('./functions/')

%% Init variables for quasiLPV
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
QuasiLPV.N = 5;
QuasiLPV.quasinorm = 1e-5;

% Simulation parameters
QuasiLPV.x0 =  [1; 0];
QuasiLPV.SimT = 20;
QuasiLPV.Tau = 0.1;

% Constraints
QuasiLPV.x_low = [-2; -5];
QuasiLPV.x_high = [2; 5];
QuasiLPV.u_low = [-1];
QuasiLPV.u_high = [1];

%Meshgrid parameters
meshgrid_res_x = 40; % Number of points in x  direction
meshgrid_res_y = 40; % Number of points in y direction

% Calculate invariant, control adm., terminal set and cost matrices
[QuasiLPV] = Calculate_sets_QuasiLPV(QuasiLPV);

%% Init variables for linear mpc

% System matrices
LinMPC.Ac = [0 1;-1 0];
LinMPC.Bc = [0; 1];
LinMPC.Cc = eye(2);
LinMPC.Dc = [0];

% Cost matrices
LinMPC.Q = [0.1 0; 0 1];
LinMPC.R = 0.01;
LinMPC.N = 5;

% Simulation parameters
LinMPC.SimT = QuasiLPV.SimT;
LinMPC.Tau = QuasiLPV.Tau;

% Constraints
LinMPC.x_low = QuasiLPV.x_low;
LinMPC.x_high = QuasiLPV.x_high;
LinMPC.u_low = QuasiLPV.u_low;
LinMPC.u_high = QuasiLPV.u_high;

% Discritize model
LinMPC.Systemc = ss(LinMPC.Ac, LinMPC.Bc, LinMPC.Cc, LinMPC.Dc);
LinMPC.Systemd = c2d(LinMPC.Systemc,LinMPC.Tau);

% Calculate invariant, control adm., terminal set and cost matrices
[LinMPC] = linMPC_Calculate_sets_LinMPC(LinMPC, []);

%% Init variables for nonlinear MPC

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

% Prepare/create nonlinear controller
NonLinMPC = perepareNonlinearMPC15(Model, MPC);

%% Prepare grid
meshgrid_xvec = linspace(QuasiLPV.x_low(1),QuasiLPV.x_high(1),meshgrid_res_x);
meshgrid_yvec = linspace(QuasiLPV.x_low(2),QuasiLPV.x_high(2),meshgrid_res_y);
statespacegrid = meshgrid(meshgrid_xvec,meshgrid_yvec);

disp('Ready to start simulating')

%% Simulate with quasiLPV

% Prepare plot
try
    close(1)
catch
    disp(' ')
end
figure(1)
hold on
grid on
xlim([QuasiLPV.x_low(1)-0.5,QuasiLPV.x_high(1)+0.5])
ylim([QuasiLPV.x_low(2)-1,QuasiLPV.x_high(2)+1])
grid on
title('Estimated feasable set quasiLPV', 'Interpreter','latex')
xlabel('State 1 (p)', 'Interpreter','latex')
ylabel('State 2 ($\dot{p}$)', 'Interpreter','latex')
plot(LinMPC.Feasiblesetx, 'Color', 'yellow', 'Alpha', 0.5);

% Simulate over all points in grid
for i = 1:length(meshgrid_xvec)
    for j = 1:length(meshgrid_yvec)
        QuasiLPV.x0 =  [meshgrid_xvec(i); meshgrid_yvec(j)];
        figure(1)
        plot(QuasiLPV.x0(1),QuasiLPV.x0(2),'*k')
        drawnow limitrate
        try
            QuasiLPV = Simulate_QuasiLPV(QuasiLPV, true);
            disp(QuasiLPV.avg_quasi_time)
        catch
            % error is renerated if infeasable
            figure(1)
            plot(QuasiLPV.x0(1),QuasiLPV.x0(2),'or')
            drawnow limitrate
        end
    end
end

% Save figure
fig1 = figure(1);
filename = fullfile('figures/', ['quasiLPV_feasable_set_estimate_' datestr(datetime,'yyyymmddHHMMSS')]);
saveas(fig1 ,filename)
saveas(fig1 ,filename, 'png')

%% Simulation of linear MPC in closed loop with Simulink model

% Prepare plot
try
    close(2)
catch
    disp(' ')
end
figure(2)
hold on
grid on
xlim([QuasiLPV.x_low(1)-0.5,QuasiLPV.x_high(1)+0.5])
ylim([QuasiLPV.x_low(2)-1,QuasiLPV.x_high(2)+1])
grid on
title('Estimated feasable set linear MPC', 'Interpreter','latex')
xlabel('State 1 (p)', 'Interpreter','latex')
ylabel('State 2 ($\dot{p}$)', 'Interpreter','latex')
plot(LinMPC.Feasiblesetx, 'Color', 'yellow', 'Alpha', 0.5);
xlim([QuasiLPV.x_low(1)-0.5,QuasiLPV.x_high(1)+0.5])
ylim([QuasiLPV.x_low(2)-1,QuasiLPV.x_high(2)+1])

% Simulate over all grid points
for i = 1:length(meshgrid_xvec)
    for j = 1:length(meshgrid_yvec)
        LinMPC.x0 =  [meshgrid_xvec(i); meshgrid_yvec(j)];
        figure(2)
        plot(LinMPC.x0(1),LinMPC.x0(2),'*k')
        drawnow limitrate
        try
            LinMPC = Simulate_NLsystem_linMPC(LinMPC, true);
            disp(LinMPC.avg_linmpc_time)
        catch
            % error is generated if infeasable
            figure(2)
            plot(LinMPC.x0(1),LinMPC.x0(2),'or')
            drawnow limitrate
        end
    end
end

% Save figure
fig2 = figure(2);
filename = fullfile('figures/', ['linMPC_feasable_set_estimate_' datestr(datetime,'yyyymmddHHMMSS')]);
saveas(fig2 ,filename)
saveas(fig2 ,filename, 'png')


%% Simulation of nonlinear MPC in closed loop with Simulink model

% Prepare plot
try
    close(3)
catch
    disp(' ')
end
figure(3)
hold on
grid on
xlim([QuasiLPV.x_low(1)-0.5,QuasiLPV.x_high(1)+0.5])
ylim([QuasiLPV.x_low(2)-1,QuasiLPV.x_high(2)+1])
grid on
title('Estimated feasable set nonlinear MPC', 'Interpreter','latex')
xlabel('State 1 (p)', 'Interpreter','latex')
ylabel('State 2 ($\dot{p}$)', 'Interpreter','latex')
plot(LinMPC.Feasiblesetx, 'Color', 'yellow', 'Alpha', 0.5);
xlim([QuasiLPV.x_low(1)-0.5,QuasiLPV.x_high(1)+0.5])
ylim([QuasiLPV.x_low(2)-1,QuasiLPV.x_high(2)+1])

% simulate over all grid points
for i = 1:length(meshgrid_xvec)
    for j = 1:length(meshgrid_yvec)
        SimSettings.x0 = [meshgrid_xvec(i); meshgrid_yvec(j)];
        figure(3)
        plot(SimSettings.x0(1),SimSettings.x0(2),'*k')
        drawnow limitrate
        try
            [NonLinMPC] = runNonlinearMPC15(SimSettings, NonLinMPC);
            disp(NonLinMPC.avg_nonlin_time)
        catch
            % Error is generated if infeasable
            figure(3)
            plot(SimSettings.x0(1),SimSettings.x0(2),'or')
            drawnow limitrate
        end
    end
end

% Save plot
fig3 = figure(3);
filename = fullfile('figures/', ['nonlinMPC_feasable_set_estimate_' datestr(datetime,'yyyymmddHHMMSS')]);
saveas(fig3 ,filename)
saveas(fig3 ,filename, 'png')

%% Print average computation times
disp(['Avg. comp. time for quasi LPV:  ',num2str(QuasiLPV.avg_quasi_time) , ' [sec]'])
disp(['Avg. comp. time for linear MPC: ',num2str(LinMPC.avg_linmpc_time) , ' [sec]'])
disp(['Avg. comp. time for nonlin MPC: ',num2str(NonLinMPC.avg_nonlin_time) , ' [sec]'])


