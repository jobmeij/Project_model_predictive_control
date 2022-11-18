
function [NonLinMPC] = runNonlinearMPC15(SimSettings,NonLinMPC)

% Set workspace to this function
Ts = SimSettings.Ts;    % Load sample time in this function workspace
x = SimSettings.x0;

% Set simulation length (samples)
NonLinMPC.Length = length(SimSettings.t);

% Pre-allocation
NonLinMPC.Utrajectory = zeros(NonLinMPC.Length,NonLinMPC.nu);
NonLinMPC.Xtrajectory = zeros(NonLinMPC.Length,NonLinMPC.nx);
NonLinMPC.usimulink = zeros(NonLinMPC.Length,NonLinMPC.nu);
NonLinMPC.x1simulink = zeros(NonLinMPC.Length,1);
NonLinMPC.x2simulink = zeros(NonLinMPC.Length,1);

% Initial conditions
NonLinMPC.x1simulink(1) = SimSettings.x0(1);
NonLinMPC.x2simulink(1) = SimSettings.x0(2);
NonLinMPC.nonlin_time = zeros(1,(NonLinMPC.Length-1));

% Run nonlinear MPC simulation
for i = 1:(NonLinMPC.Length-1)
    nonlintic = tic;
    [solutions,diagnostics] = NonLinMPC.controller{x};
    if diagnostics == 1
        error('Infeasible!');
    end
    % Input prediction horizon
    NonLinMPC.U = solutions;                                    
    NonLinMPC.Utrajectory(i) = NonLinMPC.U(1);
    
    % Use simulink simulation
    umpc = NonLinMPC.U(1);
    NonLinMPC.usimulink(i) = umpc;
    
    % Set initial conditions of integrators
    if (i == 1)
        x1_sim = SimSettings.x0(1);
        x2_sim = SimSettings.x0(2);
    else
        x1_sim = NonLinMPC.x1simulink(i);
        x2_sim = NonLinMPC.x2simulink(i);
    end
    
    % MATLAB function ODE45
%     x_sol_nl = simulateNLmodel([x1_sim; x2_sim], umpc, Ts);
%     x = x_sol_nl;
    % Simulink model
    x_sol_nl = sim_nonlinear_model([x1_sim; x2_sim], umpc, Ts);
    x = x_sol_nl';
    
    % Save results
    NonLinMPC.x1simulink(i+1) = x_sol_nl(1);
    NonLinMPC.x2simulink(i+1) = x_sol_nl(2);

    NonLinMPC.nonlin_time(i) = toc(nonlintic);
end

NonLinMPC.avg_nonlin_time = sum(NonLinMPC.nonlin_time)/(NonLinMPC.Length-1);
    
% Plot trajectories in statespace
figure(6)
plot(NonLinMPC.x1simulink,NonLinMPC.x2simulink)
grid on
title('X_1 and X_2')
xlabel('X_1')
ylabel('X_2')

xlim([NonLinMPC.x1min-0.5,NonLinMPC.x1max+0.5])
ylim([NonLinMPC.x2min-1,NonLinMPC.x2max+1])
drawnow limitrate

end
