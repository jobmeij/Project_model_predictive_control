function QuasiLPV = Simulate_NLsystem_QuasiLPV(QuasiLPV, plots)

quasinorm = 0.1;

p = sym('p', [1 QuasiLPV.N]);
% Prepare MPC QP inputs
QuasiLPV.G = (QuasiLPV.G+QuasiLPV.G')/2;

% [QuasiLPV.L_mpcsolver,~] = chol(QuasiLPV.G,'lower', 'nocheck');
% QuasiLPV.Linv = inv(QuasiLPV.L_mpcsolver);
% QuasiLPV.Linv(abs(QuasiLPV.Linv)<10*eps) = 0;

QuasiLPV.G_func = matlabFunction(QuasiLPV.G, 'vars', {p});
QuasiLPV.F_func = matlabFunction(QuasiLPV.F, 'vars', {p});
QuasiLPV.W_func = matlabFunction(QuasiLPV.W, 'vars', {p});
QuasiLPV.L_func =  matlabFunction(QuasiLPV.L, 'vars', {p});

% QuasiLPV.Linv = 
% QuasiLPV.c_p_func =  matlabFunction(QuasiLPV.G, 'vars', {x});
% QuasiLPV.b_p_func =  matlabFunction(QuasiLPV.G, 'vars', {x});
% PHI_func = matlabFunction(PHI, 'vars', {x});
% GAMMA_func = matlabFunction(GAMMA, 'vars', {x});

QuasiLPV.iA = false(size(QuasiLPV.c,1),1);
QuasiLPV.T = [0:QuasiLPV.Tau:QuasiLPV.SimT];
MPCsolverOpt = mpcqpsolverOptions();
QuasiLPV.u_sol = nan(1,length(QuasiLPV.T));
QuasiLPV.x_sol = nan(2,length(QuasiLPV.T)+1);
QuasiLPV.status = nan(1,length(QuasiLPV.T));
QuasiLPV.x_sol(:,1) = QuasiLPV.x0;

x_hor = repmat(QuasiLPV.x0,1,QuasiLPV.N);

% Simulate system
for i = 1:length(QuasiLPV.T)
    clc
    disp(['Simulating QuasiLPV - ', num2str(round(i/length(QuasiLPV.T)*100)), '%'])
    
    diff = inf;
    while(diff > quasinorm )
        
        % Solve symbolic matrices with found x
        p_sol = exp(-x_hor(1,:));
        QuasiLPV.G_sol = G_func(p_sol);
        QuasiLPV.F_sol = F_sol(p_sol);
        QuasiLPV.W_sol = W_sol(p_sol);
        QuasiLPV.L_sol = L_sol(p_sol);
        
        [QuasiLPV.L_mpcsolver,~] = chol(QuasiLPV.G_func,'lower', 'nocheck');
        QuasiLPV.Linv = inv(QuasiLPV.L_mpcsolver);
        QuasiLPV.Linv(abs(QuasiLPV.Linv)<10*eps) = 0;
        
        
        QuasiLPV.c_p = QuasiLPV.F_func*QuasiLPV.x_sol(:,i);
        QuasiLPV.b_p = QuasiLPV.c+QuasiLPV.W_func*QuasiLPV.x_sol(:,i);
        [U_solve,QuasiLPV.status(i),QuasiLPV.iA,~] = mpcqpsolver(QuasiLPV.Linv,QuasiLPV.c_p_func, -QuasiLPV.L_func, -QuasiLPV.b_p, [], zeros(0,1) ,QuasiLPV.iA, MPCsolverOpt);
        if QuasiLPV.status(i) < 0 % unfeasable
            QuasiLPV.feasable = false;
            error('Infeasable solution')
        end
        U_solve_prev = U_solve;
        for j = 1:N
            % MATLAB function ODE45
            QuasiLPV.x_hor(:,j+1) = simulateNLmodel(QuasiLPV.x_hor(:,j), U_solve(j), QuasiLPV.Tau);
            % Simulink model
            % QuasiLPV.x_hor(:,j+1) = sim_nonlinear_model(QuasiLPV.x_hor(:,j), U_solve(j), QuasiLPV.Tau);
        end
        diff = abs((U_solve_prev - U_solve).^2)
    end
    
    QuasiLPV.u_sol(i) = U_solve(1,:);
    QuasiLPV.x_sol(:,i+1) = QuasiLPV.x_hor(:,2);
    
    % Prepare x horizon for next iteration by shifting it one
    x_hor = [x_hor(:,2:end), x_hor(:,end)];
end
if i == length(QuasiLPV.T_nl)
    QuasiLPV.feasable = true;
end

if plots == true
    % Plot results
    figure(4)
    subplot(2,1,1)
    plot(QuasiLPV.T_nl, QuasiLPV.x_sol_nl(1,1:(end-1)), QuasiLPV.T_nl, QuasiLPV.x_sol_nl(2,1:(end-1)))
    grid on
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [sec]', 'Interpreter','latex');
    title('States over time')
    legend('State 1 (p)', 'State 2 ($\dot{p}$)', 'location', 'best', 'Interpreter','latex')
    subplot(2,1,2)
    plot(QuasiLPV.T_nl,QuasiLPV.u_sol_nl(1:end))
    grid on
    title('Control input u(t)')
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [sec]', 'Interpreter','latex');
    sgtitle('Nonlinear continious time simulation', 'Interpreter','latex');

    figure(5)
    plot(QuasiLPV.Feasiblesetx, 'Color', 'Yellow'); 
    hold on
    plot(QuasiLPV.InvSetXU,'Color', 'Green', 'Alpha', 0.5); 
    plot(QuasiLPV.x_sol_nl(1,:), QuasiLPV.x_sol_nl(2,:), '-Blue')
    plot(QuasiLPV.x_sol_nl(1,1), QuasiLPV.x_sol_nl(2,1), 'ored')
    plot(QuasiLPV.x_sol_nl(1,i), QuasiLPV.x_sol_nl(2,i), '+red')
    title('Simulation result of the nonlinear continious time plant', 'Interpreter','latex')
    xlabel('State 1 (p)', 'Interpreter','latex')
    ylabel('State 2 ($\dot{p}$)', 'Interpreter','latex')
    legend('Feasable set (N=50)','Invariant set', 'Solution', 'Init condition', 'Final state', 'location','Northwest')
end