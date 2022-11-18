function LinMPC = Simulate_NLsystem_linMPC(LinMPC, plots)

% Prepare MPC QP inputs
LinMPC.G=(LinMPC.G+LinMPC.G')/2;
[LinMPC.L_mpcsolver,~] = chol(LinMPC.G,'lower');
LinMPC.Linv = inv(LinMPC.L_mpcsolver);
LinMPC.Linv(abs(LinMPC.Linv)<10*eps) = 0;
LinMPC.iA = false(size(LinMPC.c,1),1);
LinMPC.T_nl = [0:LinMPC.Tau:LinMPC.SimT];
MPCsolverOpt = mpcqpsolverOptions();
LinMPC.u_sol_nl = nan(1,length(LinMPC.T_nl));
LinMPC.x_sol_nl = nan(2,length(LinMPC.T_nl)+1);
LinMPC.status_nl = nan(1,length(LinMPC.T_nl));
LinMPC.x_sol_nl(:,1) = LinMPC.x0;
LinMPC.linmpc_time = zeros(1,length(LinMPC.T_nl));

% Simulate system
for i = 1:length(LinMPC.T_nl)
    lintic = tic;
    clc
    disp(['Simulating Nonlinear model - ', num2str(round(i/length(LinMPC.T_nl)*100)), '%'])
    LinMPC.c_p = LinMPC.F*LinMPC.x_sol_nl(:,i);
    LinMPC.b_p = LinMPC.c+LinMPC.W*LinMPC.x_sol_nl(:,i);
    [U_solve,LinMPC.status_nl(i),LinMPC.iA,~] = mpcqpsolver(LinMPC.Linv,LinMPC.c_p, -LinMPC.L, -LinMPC.b_p, [], zeros(0,1) ,LinMPC.iA, MPCsolverOpt);
    if LinMPC.status_nl(i) < 0 % unfeasable
            LinMPC.feasable = false;
            error('Infeasable solution')
    end
    if i == 1
        LinMPC.U_Init_hor = U_solve;
    end
    LinMPC.u_sol_nl(i) = U_solve(1,:);
    
    % MATLAB function ODE45
%     LinMPC.x_sol_nl(:,i+1) = simulateNLmodel(LinMPC.x_sol_nl(:,i), LinMPC.u_sol_nl(i), LinMPC.Tau);
    % Simulink model
    LinMPC.x_sol_nl(:,i+1) = sim_nonlinear_model(LinMPC.x_sol_nl(:,i), LinMPC.u_sol_nl(i), LinMPC.Tau);
    
    LinMPC.linmpc_time(i) = toc(lintic);
end
LinMPC.feasable = true;

% Calculate average computation time
LinMPC.avg_linmpc_time = sum(LinMPC.linmpc_time)/length(LinMPC.T_nl);

if plots == true
    % Plot results
    figure(5)
    plot(LinMPC.x_sol_nl(1,:), LinMPC.x_sol_nl(2,:), '-Blue')
    hold on
    plot(LinMPC.x_sol_nl(1,1), LinMPC.x_sol_nl(2,1), 'ored')
    plot(LinMPC.x_sol_nl(1,i), LinMPC.x_sol_nl(2,i), '+red')
    hold off
    grid on

    xlabel('State 1 (p)', 'Interpreter','latex')
    ylabel('State 2 ($\dot{p}$)', 'Interpreter','latex')
    xlim([LinMPC.x_low(1)-0.5,LinMPC.x_high(1)+0.5])
    ylim([LinMPC.x_low(2)-1,LinMPC.x_high(2)+1])
    drawnow limitrate
end