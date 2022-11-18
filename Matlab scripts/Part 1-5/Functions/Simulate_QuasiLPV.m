function QuasiLPV = Simulate_QuasiLPV(QuasiLPV, plots)

% Prepare MPC QP inputs
QuasiLPV.iA = false(size(QuasiLPV.c,1),1);
QuasiLPV.T = [0:QuasiLPV.Tau:QuasiLPV.SimT];
MPCsolverOpt = mpcqpsolverOptions();
QuasiLPV.u_sol = nan(1,length(QuasiLPV.T));
QuasiLPV.x_sol = nan(2,length(QuasiLPV.T)+1);
QuasiLPV.status = nan(1,length(QuasiLPV.T));
QuasiLPV.x_sol(:,1) = QuasiLPV.x0;

x_hor = repmat(QuasiLPV.x0,1,QuasiLPV.N);
QuasiLPV.feasable = true;
QuasiLPV.quasi_time = zeros(1,length(QuasiLPV.T));

% Simulate system
for i = 1:length(QuasiLPV.T)
    tic_quasi = tic;
    clc
    disp(['Simulating QuasiLPV - ', num2str(round(i/length(QuasiLPV.T)*100)), '%'])
    
    diff = inf;
    U_solve_prev = inf(QuasiLPV.N,1);
    while(diff > QuasiLPV.quasinorm )
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
            error('Infeasable solution')
        end
       
        for j = 1:QuasiLPV.N
            % MATLAB function ODE45
%             x_hor(:,j+1) = simulateNLmodel(x_hor(:,j), U_solve(j), QuasiLPV.Tau);
            % Simulink model
            x_hor(:,j+1) = sim_nonlinear_model(x_hor(:,j), U_solve(j), QuasiLPV.Tau);
        end
        % Calculate difference for quasi norm
        diff = norm(U_solve_prev - U_solve);
        U_solve_prev = U_solve;
    end
    
    QuasiLPV.u_sol(i) = U_solve(1,:);
    QuasiLPV.x_sol(:,i+1) = x_hor(:,2);
    
    % Prepare x horizon for next iteration by shifting it one
    x_hor = [x_hor(:,2:end), x_hor(:,end)];
    QuasiLPV.quasi_time(i) = toc(tic_quasi);
end

% Calculate avg. computation time
QuasiLPV.avg_quasi_time = sum(QuasiLPV.quasi_time)/length(QuasiLPV.T);


if plots == true
    % Plot results
    figure(4)
    plot(QuasiLPV.x_sol(1,:), QuasiLPV.x_sol(2,:), '-Blue')
    xlim([QuasiLPV.x_low(1)-1,QuasiLPV.x_high(1)+1])
    ylim([QuasiLPV.x_low(2)-1,QuasiLPV.x_high(2)+1])
    hold on
    plot(QuasiLPV.x_sol(1,1), QuasiLPV.x_sol(2,1), 'ored')
    plot(QuasiLPV.x_sol(1,i), QuasiLPV.x_sol(2,i), '+red')
    hold off
    drawnow limitrate
end