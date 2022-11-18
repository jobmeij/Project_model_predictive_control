% MPC graded homework assignment
% Part 1.1 -  initial nonlinear system simulation
% Job Meijer            - 1268155
% Marcel van Wensveen   - 1253085

clear all; close all; clc;

% Add functions directory
addpath('./functions/')

% Simulation settings
Ts = 1/10;                      % Simulation stop time [s]
x0 = [1; 0];                    % Initial conditions
Tstep = Ts;                     % Step time of simulation [s]
Tend = 20;                      % End time of simulation [s]
t = 0:Tstep:Tend;               % Simulation time vector
x_sol = zeros(2,length(t));     % Preallocate memory
x_sol(:,1) = x0;                % Start with initial conditions


%% Run simulation model

 for k = 1:(length(t)-1)
    clc
    disp(['Simulating Nonlinear model - ', num2str(round(k/length(t)*100)), '%'])
    x_sol(:,k+1) = sim_nonlinear_model(x_sol(:,k), 0, Ts);
 end

%% Plotting results
figure(1)
set(gcf,'Position',[500 300 1200 500])

% 1D plot
subplot (2,3,1)         
hold on
plot(t,x_sol(1,:))
plot(t,x_sol(2,:))
hold off
grid on
xlabel('Time [s]', 'Interpreter','latex')
ylabel('Amplitude [-]', 'Interpreter','latex')
legend('${p}(t)$','$\dot{p}(t)$','Location','Best','Interpreter','latex')
title('1D plot of simulation over time', 'Interpreter','latex')

% 2D plot
subplot (2,3,4)         
plot(x_sol(1,:),x_sol(2,:))
grid on
title('Statespace plot of simulation', 'Interpreter','latex')
xlabel('${p}(t)$', 'Interpreter','latex')
ylabel('$\dot{p}(t)$', 'Interpreter','latex')

% 3D plot
subplot (2,3,[2,3,5,6])     
plot3(t,x_sol(1,:),x_sol(2,:))
set(gca,'ydir','reverse')
view(45,-10)
grid on
title('3D plot of simulation', 'Interpreter','latex')
xlabel('Time [s]', 'Interpreter','latex')
ylabel('${p}(t)$', 'Interpreter','latex')
zlabel('$\dot{p}(t)$', 'Interpreter','latex')

