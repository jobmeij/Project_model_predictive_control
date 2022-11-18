% MPC graded homework assignment
% Part 1.2 -  MPC design for constant set-point regulation
% Job Meijer            - 1268155
% Marcel van Wensveen   - 1253085


% Init script
clear all
close all
close_system    % Close all open simulink models in memory
clc

% Add functions directory
addpath('./functions/')

%% Init variables for simulation

% Uncomment the desired conditions for the simulation
% linearizationpoint = '00';
linearizationpoint = '10';
% costset = 'GivenCost';
costset = 'NonlinearCost';


% pich the selected conditions
switch linearizationpoint
    case  '10'
        LinMPC.Ac = [0 1;-1 0];
        LinMPC.Bc = [0; 1];
        LinMPC.Cc = eye(2);
        LinMPC.Dc = [0];
    case '00'
        LinMPC.Ac = [0 1;-1 2];
        LinMPC.Bc = [0; 1];
        LinMPC.Cc = eye(2);
        LinMPC.Dc = [0];
end

switch costset
    case 'GivenCost'
        LinMPC.Q = [1 0; 0 4];
        LinMPC.R = 0.1;
        LinMPC.N = 50;
    case 'NonlinearCost'
        LinMPC.Q = [0.1 0; 0 1];
        LinMPC.R = 0.01;
        LinMPC.N = 50;
end

% Simulation parameters
LinMPC.x0 =  [1; 0];
LinMPC.SimT = 20;
LinMPC.Tau = 0.1;

% Constraints
LinMPC.x_low = [-2; -5];
LinMPC.x_high = [2; 5];
LinMPC.u_low = [-1];
LinMPC.u_high = [1];

% Discritize model
LinMPC.Systemc = ss(LinMPC.Ac, LinMPC.Bc, LinMPC.Cc, LinMPC.Dc);
LinMPC.Systemd = c2d(LinMPC.Systemc,LinMPC.Tau);

%% Calculate invariant, control adm., terminal set and cost matrices
[LinMPC] = Calculate_sets_LinMPC_12(LinMPC, 1);

%% Simulate with linear discrete time model
LinMPC = Simulate_linear_MPC_12(LinMPC, true);

%% Simulate with nonlinear continious time model
LinMPC = Simulate_NLsystem_linMPC_12(LinMPC, true);

%% plots for report

figure(7)
subplot(2,1,1)
plot(LinMPC.T_lin, LinMPC.x_sol_lin(1,1:(end-1)), LinMPC.T_lin, LinMPC.x_sol_lin(2,1:(end-1)))
hold on
plot(LinMPC.T_nl, LinMPC.x_sol_nl(1,1:(end-1)), LinMPC.T_nl, LinMPC.x_sol_nl(2,1:(end-1)))
hold off
grid on
ylabel('Amplitude', 'Interpreter','latex');
xlabel('Time [sec]', 'Interpreter','latex');
title('States over time')
legend('State 1 Lin system', 'State 2 Lin system','State 1 NL system', 'State 2 NL system', 'Interpreter','latex')
subplot(2,1,2)
plot(LinMPC.T_lin,LinMPC.u_sol_lin(1:end))
hold on
plot(LinMPC.T_nl,LinMPC.u_sol_nl(1:end))
hold off
grid on
title('Control input u(t)')
legend('Lin system', 'NL system', 'Interpreter','latex')
ylabel('Amplitude', 'Interpreter','latex');
xlabel('Time [sec]', 'Interpreter','latex');
sgtitle('Linear MPC simulation', 'Interpreter','latex');

figure(8)
plot(LinMPC.Feasiblesetx, 'Color', 'Yellow');
hold on
plot(LinMPC.InvSetXU,'Color', 'Green', 'Alpha', 0.5);
plot(LinMPC.x_sol_lin(1,:), LinMPC.x_sol_lin(2,:), '-Blue')
plot(LinMPC.x_sol_nl(1,:), LinMPC.x_sol_nl(2,:), '-red')
title('Simulation result of the linear MPC', 'Interpreter','latex')
xlabel('State 1 (p)', 'Interpreter','latex')
ylabel('State 2 ($\dot{p}$)', 'Interpreter','latex')
legend('Feasable set (N=50)','Invariant set', 'Solution Lin system', 'Solution NL system', 'location','Northwest')


