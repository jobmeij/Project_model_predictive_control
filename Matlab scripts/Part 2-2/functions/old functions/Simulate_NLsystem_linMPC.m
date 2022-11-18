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

% Simulate system
for i = 1:length(LinMPC.T_nl)
    clc
    disp(['Simulating Nonlinear model - ', num2str(round(i/length(LinMPC.T_nl)*100)), '%'])
    LinMPC.c_p = LinMPC.F*LinMPC.x_sol_nl(:,i);
    LinMPC.b_p = LinMPC.c+LinMPC.W*LinMPC.x_sol_nl(:,i);
    [U_solve,LinMPC.status_nl(i),LinMPC.iA,~] = mpcqpsolver(LinMPC.Linv,LinMPC.c_p, -LinMPC.L, -LinMPC.b_p, [], zeros(0,1) ,LinMPC.iA, MPCsolverOpt);
    if LinMPC.status_nl(i) < 0 % unfeasable
            LinMPC.feasable = false;
            disp('Infeasable solution')
            break
    end
    if i == 1
        LinMPC.U_Init_hor = U_solve;
    end
    LinMPC.u_sol_nl(i) = U_solve(1,:);
    LinMPC.x_sol_nl(:,i+1) = sim_nonlinear_model(LinMPC.x_sol_nl(:,i), LinMPC.u_sol_nl(i), LinMPC.Tau);
end
if i == length(LinMPC.T_nl)
    LinMPC.feasable = true;
end

if plots == true
    % Plot results
    figure(4)
    subplot(2,1,1)
    plot(LinMPC.T_nl, LinMPC.x_sol_nl(1,1:(end-1)), LinMPC.T_nl, LinMPC.x_sol_nl(2,1:(end-1)))
    grid on
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [sec]', 'Interpreter','latex');
    title('States over time')
    legend('State 1 (p)', 'State 2 ($\dot{p}$)', 'location', 'best', 'Interpreter','latex')
    subplot(2,1,2)
    plot(LinMPC.T_nl,LinMPC.u_sol_nl(1:end))
    grid on
    title('Control input u(t)')
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [sec]', 'Interpreter','latex');
    sgtitle('Nonlinear continious time simulation', 'Interpreter','latex');

    figure(5)
    plot(LinMPC.Feasiblesetx, 'Color', 'Yellow'); 
    hold on
    plot(LinMPC.InvSetXU,'Color', 'Green', 'Alpha', 0.5); 
    plot(LinMPC.x_sol_nl(1,:), LinMPC.x_sol_nl(2,:), '-Blue')
    plot(LinMPC.x_sol_nl(1,1), LinMPC.x_sol_nl(2,1), 'ored')
    plot(LinMPC.x_sol_nl(1,i), LinMPC.x_sol_nl(2,i), '+red')
    title('Simulation result of the nonlinear continious time plant', 'Interpreter','latex')
    xlabel('State 1 (p)', 'Interpreter','latex')
    ylabel('State 2 ($\dot{p}$)', 'Interpreter','latex')
    legend('Feasable set (N=50)','Invariant set', 'Solution', 'Init condition', 'Final state', 'location','Northwest')
end