%% Part 1.3: Closed loop simulation with linear MPC controller
function [LinMPC] = runLinearMPC(SimSettings,Model,MPC)
x = SimSettings.x0';                  % Initial state conditions

% Set workspace to this function
options = simset('SrcWorkspace','current');
Ts = SimSettings.Ts;    % Load sample time in this function workspace

% Create input of zeros
LinMPC.u_sim = zeros(length(SimSettings.t),1);

% Pre-allocate
LinMPC.x1_sim = zeros(length(SimSettings.t),1);
LinMPC.x2_sim = zeros(length(SimSettings.t),1);

% Initial conditions
LinMPC.x1_sim = SimSettings.x0(1);
LinMPC.x2_sim = SimSettings.x0(2);

% Controller options
LinMPC.EnableQp = true;
LinMPC.EnableInputConstraints = false; %true
LinMPC.EnableStateConstraints = false; %true
LinMPC.options_qp =  optimoptions('quadprog','Display','off');     % quadprog option

% Set simulation length
LinMPC.Length = length(SimSettings.t);

% Run simulation loop
for k = 1:(LinMPC.Length-1)
    % Set input
    umpc = LinMPC.u_sim(k);        % Use given input for simulation
    
    % Overwrite input with QP input
    if LinMPC.EnableQp
        [LinMPC.u_qp,LinMPC.fval,LinMPC.exitflag] = quadprog(MPC.G,MPC.F*x',MPC.L,MPC.c+MPC.W*x',[],[],[],[],[],LinMPC.options_qp);
%         if LinMPC.exitflag ~= 1
%             warning('exitflag quadprog =%d\n', LinMPC.exitflag)
            if LinMPC.exitflag == -2
                error('Optimization problem is infeasible.')
            end
%         end
        umpc = LinMPC.u_qp(1);
    end
    
    % Input constraints (if enabled)
    if LinMPC.EnableInputConstraints
        if (umpc < Model.uMin)
            umpc = Model.uMin;
        elseif (umpc > Model.uMax)
            umpc = Model.uMax;
        end
    end
    
    % Run simulation
    sim('CTmodel.slx',[],options);
    
    % Save simulation results for plotting
    LinMPC.t_sim(k) = Ts*(k-1);
    LinMPC.u_sim(k) = umpc;
    LinMPC.x1_sim(k+1) = x1simu(end);    % Select end value
    LinMPC.x2_sim(k+1) = x2simu(end);    % Select end value

    % State constraints (if enabled)
    if LinMPC.EnableStateConstraints
        if (x1simu(end) < Model.xMin(1))
            x1simu(end) = Model.xMin(1);
        elseif (x1simu(end) > Model.xMax(1))
            x1simu(end) = Model.xMax(1);
        end
        
        if (x2simu(end) < Model.xMin(2))
            x2simu(end) = Model.xMin(2);
        elseif (x2simu(end) > Model.xMax(2))
            x2simu(end) = Model.xMax(2);
        end
    end
    
    % Set initial conditions for next iteration
    x = [LinMPC.x1_sim(k+1) LinMPC.x2_sim(k+1)];
end

if (SimSettings.LinMPCplot)
    % 1D plotting
    figure(5)
    set(gcf,'Position',[10 100 1200 400])
    hold on
    plot(LinMPC.t_sim,LinMPC.u_sim)
    plot(LinMPC.t_sim,LinMPC.x1_sim)
    plot(LinMPC.t_sim,LinMPC.x2_sim)
    hold off
    grid on
    xlabel('Time [s]', 'Interpreter','latex')
    ylabel('Amplitude [-]', 'Interpreter','latex')
    legend('$u_{mpc}$','$\dot{p}(t)$','$\ddot{p}(t)$','Location','Best','Interpreter','latex')
    title('Continuous-time simulation', 'Interpreter','latex')
end

end