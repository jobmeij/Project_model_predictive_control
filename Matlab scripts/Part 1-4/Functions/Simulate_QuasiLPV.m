function QuasiLPV = Simulate_QuasiLPV(QuasiLPV, plots)

quasinorm = 1e-5;

p = sym('p', [1 QuasiLPV.N]);
% Prepare MPC QP inputs
QuasiLPV.G = (QuasiLPV.G+QuasiLPV.G')/2;

QuasiLPV.G_func = matlabFunction(QuasiLPV.G, 'vars', {p});
QuasiLPV.F_func = matlabFunction(QuasiLPV.F, 'vars', {p});
QuasiLPV.W_func = matlabFunction(QuasiLPV.W, 'vars', {p});
QuasiLPV.L_func =  matlabFunction(QuasiLPV.L, 'vars', {p});

QuasiLPV.iA = false(size(QuasiLPV.c,1),1);
QuasiLPV.T = [0:QuasiLPV.Tau:QuasiLPV.SimT];
MPCsolverOpt = mpcqpsolverOptions();
QuasiLPV.u_sol = nan(1,length(QuasiLPV.T));
QuasiLPV.x_sol = nan(2,length(QuasiLPV.T)+1);
QuasiLPV.status = nan(1,length(QuasiLPV.T));
QuasiLPV.x_sol(:,1) = QuasiLPV.x0;

x_hor = repmat(QuasiLPV.x0,1,QuasiLPV.N);
QuasiLPV.feasable = true;
% Simulate system
for i = 1:length(QuasiLPV.T)
    clc
    disp(['Simulating QuasiLPV - ', num2str(round(i/length(QuasiLPV.T)*100)), '%'])
    
    diff = inf;
    U_solve_prev = inf(QuasiLPV.N,1);
    while(diff > quasinorm )
        
        % Solve symbolic matrices with found x
        p_sol = QuasiLPV.p_eq(x_hor(1,:));
        QuasiLPV.G_sol = QuasiLPV.G_func(p_sol);
        QuasiLPV.F_sol = QuasiLPV.F_func(p_sol);
        QuasiLPV.W_sol = QuasiLPV.W_func(p_sol);
        QuasiLPV.L_sol = QuasiLPV.L_func(p_sol);
        
        [QuasiLPV.L_mpcsolver,~] = chol(QuasiLPV.G_sol,'lower');
        QuasiLPV.Linv = inv(QuasiLPV.L_mpcsolver);
        QuasiLPV.Linv(abs(QuasiLPV.Linv)<10*eps) = 0;
        
        
        QuasiLPV.c_p = QuasiLPV.F_sol*QuasiLPV.x_sol(:,i);
        QuasiLPV.b_p = QuasiLPV.c+QuasiLPV.W_sol*QuasiLPV.x_sol(:,i);
        [U_solve,QuasiLPV.status(i),QuasiLPV.iA,~] = mpcqpsolver(QuasiLPV.Linv,QuasiLPV.c_p, -QuasiLPV.L_sol, -QuasiLPV.b_p, [], zeros(0,1) ,QuasiLPV.iA, MPCsolverOpt);
        if QuasiLPV.status(i) < 0 % unfeasable
            QuasiLPV.feasable = false;
            disp('Infeasable solution')
            U_solve = max(min(U_solve,QuasiLPV.u_high),QuasiLPV.u_low);
            error('unfeasable')
        end
       
        for j = 1:QuasiLPV.N
            % MATLAB function ODE45
%             x_hor(:,j+1) = simulateNLmodel(x_hor(:,j), U_solve(j), QuasiLPV.Tau);
            % Simulink model
            x_hor(:,j+1) = sim_nonlinear_model(x_hor(:,j), U_solve(j), QuasiLPV.Tau);
        end
        diff = norm(U_solve_prev - U_solve);
        U_solve_prev = U_solve;
    end
    
    QuasiLPV.u_sol(i) = U_solve(1,:);
    QuasiLPV.x_sol(:,i+1) = x_hor(:,2);
    
    % Prepare x horizon for next iteration by shifting it one
    x_hor = [x_hor(:,2:end), x_hor(:,end)];
end

if plots == true
    % Plot results
    figure(4)
    subplot(2,1,1)
    plot(QuasiLPV.T, QuasiLPV.x_sol(1,1:(end-1)), QuasiLPV.T, QuasiLPV.x_sol(2,1:(end-1)))
    grid on
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [sec]', 'Interpreter','latex');
    title('States over time')
    legend('State 1 (p)', 'State 2 ($\dot{p}$)', 'location', 'best', 'Interpreter','latex')
    subplot(2,1,2)
    plot(QuasiLPV.T,QuasiLPV.u_sol(1:end))
    grid on
    title('Control input u(t)')
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [sec]', 'Interpreter','latex');
    sgtitle('QuasiLPV simulation', 'Interpreter','latex');

    figure(5)
    plot(QuasiLPV.x_sol(1,:), QuasiLPV.x_sol(2,:), '-Blue')
    hold on
    plot(QuasiLPV.x_sol(1,1), QuasiLPV.x_sol(2,1), 'ored')
    plot(QuasiLPV.x_sol(1,i), QuasiLPV.x_sol(2,i), '+red')
    title('QuasiLPV simulation in statespace', 'Interpreter','latex')
    xlabel('State 1 (p)', 'Interpreter','latex')
    ylabel('State 2 ($\dot{p}$)', 'Interpreter','latex')
    legend('Solution', 'Init condition', 'Final state', 'location','best')
end