function refMPC = Simulate_refMPC(refMPC, plots)

% Prepare MPC QP inputs
refMPC.G=(refMPC.G+refMPC.G')/2;
[refMPC.L_mpcsolver,~] = chol(refMPC.G,'lower');
refMPC.Linv = inv(refMPC.L_mpcsolver);
refMPC.Linv(abs(refMPC.Linv)<100*eps) = 0;
refMPC.iA = false(size(refMPC.c,1),1);
refMPC.T_lin = [0:refMPC.Tau:refMPC.SimT];
MPCsolverOpt = mpcqpsolverOptions();
refMPC.u_sol_lin = nan(1,length(refMPC.T_lin));
refMPC.x_sol_lin = nan(2,length(refMPC.T_lin)+1);

refMPC.status_lin = nan(1,length(refMPC.T_lin));
refMPC.x_sol_lin(:,1) = refMPC.x0;

% Simulate system
for i = 1:length(refMPC.T_lin)
    % Display progress
    clc
    disp(['Simulating linear model - ', num2str(round(i/length(refMPC.T_lin)*100)), '%'])
    
    % Estimated state values with disturbance removed
    refMPC.c_p = refMPC.F*(refMPC.x_dist_est(:,i)-refMPC.x_ss(:,i)) + refMPC.Y*refMPC.ref(i:(i+refMPC.N));
    refMPC.b_p = refMPC.c + refMPC.W*(refMPC.x_dist_est(:,i)-refMPC.x_ss(:,i));
    
    % Solve MPC problem
    [U_solve,refMPC.status_lin(i),refMPC.iA,~] = mpcqpsolver(refMPC.Linv, refMPC.c_p, -refMPC.L, -refMPC.b_p, [], zeros(0,1) ,refMPC.iA, MPCsolverOpt);
    if refMPC.status_lin(i) < 0 % unfeasable
            refMPC.feasable = false;
            disp('Infeasable solution')
            break
    end
    
    % Initial condition
    if i == 1
        refMPC.U_Init_hor = U_solve;
    end
    
    % Remove steady state u from input
    refMPC.u_sol_lin(i) = (U_solve(1,:) + refMPC.u_ss(i));
    
    % Calculate states and output
    refMPC.x_sol_lin(:,i+1) = refMPC.Ad*refMPC.x_sol_lin(:,i) + refMPC.Bd*refMPC.u_sol_lin(i) + refMPC.BdDist*refMPC.dk(i);
    refMPC.y_sol_lin(i+1) = refMPC.Cd*refMPC.x_sol_lin(:,i+1);

    % Disturbance estimation    
    refMPC.x_dist_est(:,i+1) = (refMPC.Ad-refMPC.Lx_dist*refMPC.Cd)*refMPC.x_dist_est(:,i) + refMPC.BdDist*refMPC.d_dist_est(i) + refMPC.Bd*refMPC.u_sol_lin(i) + refMPC.Lx_dist*refMPC.y_sol_lin(i);   
    refMPC.d_dist_est(i+1) = -refMPC.Ld_dist*refMPC.Cd*refMPC.x_dist_est(:,i) + refMPC.d_dist_est(i) + refMPC.Ld_dist*refMPC.y_sol_lin(i);
    
    % Computing x steady state and u steady state
    XU_ss = inv([refMPC.Ad-eye(2), refMPC.Bd; refMPC.Cd, 0])*[-refMPC.BdDist*refMPC.d_dist_est(i); 0];
    refMPC.x_ss(1,i+1) = XU_ss(1);
    refMPC.x_ss(2,i+1) = XU_ss(2);
    refMPC.u_ss(i+1) = XU_ss(3);

end

% Check is simulation finished succesfully
if i == length(refMPC.T_lin)
    refMPC.feasable = true;
end

% Plotting
if plots
    % Plot results
    figure(1)
    subplot(3,1,1)
    plot(refMPC.T_lin, refMPC.x_sol_lin(1,1:(end-1)), refMPC.T_lin, refMPC.ref(1:(end-refMPC.N)))
    grid on
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [k]', 'Interpreter','latex');
    title('Angle of elevation (\alpha)')
    legend('Actual position', 'Reference', 'location', 'best', 'Interpreter','latex')
    subplot(3,1,2)
    plot(refMPC.T_lin, refMPC.x_sol_lin(2,1:(end-1)))
    grid on
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [k]', 'Interpreter','latex');
    title('Pitch rate (q)')
    subplot(3,1,3)
    plot(refMPC.T_lin,refMPC.u_sol_lin(1:end))
    grid on
    title('Control input u(k)')
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [k]', 'Interpreter','latex');
    sgtitle('Longitudinal model of an aircraft', 'Interpreter','latex');
       
    figure(2)
    hold on
    plot(refMPC.T_lin, refMPC.x_ss(1,1:(end-1)))
    plot(refMPC.T_lin, refMPC.x_ss(2,1:(end-1)))
    plot(refMPC.T_lin, refMPC.u_ss(1,1:(end-1)))
    hold off
    grid on
    xlabel('Time [k]')
    legend('x^s^s 1','x^s^s 2','u^s^s')
    title('x^s^s and u^s^s')

    figure(3)
    hold on
    % states
    plot(refMPC.T_lin, refMPC.x_dist_est(1,1:(end-1)))
    plot(refMPC.T_lin, refMPC.x_dist_est(2,1:(end-1)))
    plot(refMPC.T_lin, refMPC.x_sol_lin(1,1:(end-1)),'--')
    plot(refMPC.T_lin, refMPC.x_sol_lin(2,1:(end-1)),'--')
    % disturbance
    plot(refMPC.T_lin, refMPC.d_dist_est(1,1:(end-1)),'LineWidth',2)
    plot(refMPC.T_lin, refMPC.dk(1:(end),1),'--','LineWidth',2)
    hold off
    grid on
    xlabel('Time [k]')
    title('Inserted and estimated disturbance')
    legend('x_e_s_t1','x_e_s_t2','x1','x2','d_e_s_t','d')
    
    % Plot for report:
    figure(4)  
    set(gcf,'position',[10,100,1600,600])
    subplot 121
    hold on
    plot(refMPC.T_lin, refMPC.ref(1:(end-refMPC.N),1),'LineWidth',2,'Color','black')
    plot(refMPC.T_lin, refMPC.x_sol_lin(1,1:(end-1)),'--','LineWidth',2,'Color','red')  
    hold off
    grid on
    ylabel('Amplitude', 'Interpreter','latex');
    xlabel('Time [k]', 'Interpreter','latex');
    title('Reference and measured angle of elevation','Interpreter','latex')
    legend('Reference', 'Measured', 'location', 'best', 'Interpreter','latex')
    
    subplot 122
    hold on
    plot(refMPC.T_lin, refMPC.dk(1:(end),1),'LineWidth',2,'Color','black')
    plot(refMPC.T_lin, refMPC.d_dist_est(1,1:(end-1)),'--','LineWidth',2,'Color','red')
    hold off
    grid on
    xlabel('Time [k]','Interpreter','latex')
    title('Inserted and estimated disturbance d(k)','Interpreter','latex')
    legend('Inserted','Estimated','Interpreter','latex')
end