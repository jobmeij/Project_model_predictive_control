%% Part 1.3: Simulating the CT model with given initial condition and zero input
function [CTsim, MPC] = runContinuousTimeModel(SimSettings, Model, MPC)
x = SimSettings.x0';                  % Initial state conditions

% Set workspace to this function
options = simset('SrcWorkspace','current');
Ts = SimSettings.Ts;    % Load sample time in this function workspace

% Create input
CTsim.u_sim = zeros(length(SimSettings.t),1);

% Pre-allocate
CTsim.x1_sim = zeros(length(SimSettings.t),1);
CTsim.x2_sim = zeros(length(SimSettings.t),1);

% Controller options
CTsettings.EnableQp = false;
CTsettings.EnableInputConstraints = false;
CTsettings.EnableStateConstraints = false;
CTsettings.options_qp =  optimoptions('quadprog','Display','off');     % quadprog option

% Run simulation loop
for k = 1:length(SimSettings.t)
    % Set input
    umpc = CTsim.u_sim(k);        % Use given input for simulation

    % Overwrite input with QP input
    if CTsettings.EnableQp
        [CTsim.u_qp,CTsim.fval,CTsim.exitflag] = quadprog(MPC.G,MPC.F*CTsim.x',MPC.L,MPC.c+MPC.W*CTsim.x',[],[],[],[],[],CTsettings.options_qp);
        if CTsim.exitflag ~= 1
            warning('exitflag quadprog =%d\n', CTsim.exitflag)
            if CTsim.exitflag == -2
                fprint('Optimization problem is infeasible. \n')
            end
        end
        umpc = CTsim.u_qp(1);
    end
    
    % Input constraints (if enabled)
    if CTsettings.EnableInputConstraints
        if (umpc < Model.uMin)
            umpc = Model.uMin;
        elseif (umpc > Model.uMax)
            umpc = Model.uMax;
        end
    end
%     umpc
    
    % Run Simulink simulation
    sim('CTmodel.slx',[],options);

    % State constraints (if enabled)
    if CTsettings.EnableStateConstraints
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
    
    % Save simulation results for plotting
    CTsim.t_sim(k) = Ts*(k-1);
    CTsim.u_sim(k) = umpc;
    CTsim.x1_sim(k) = x1simu(end);    % Select end value
    CTsim.x2_sim(k) = x2simu(end);    % Select end value
    
    % Set initial conditions for next iteration
    x = [CTsim.x1_sim(k) CTsim.x2_sim(k)];
end

if (SimSettings.CTplot)
    figure()
    hold on
    plot(1:length(SimSettings.t),CTsim.x1_sim)
    plot(1:length(SimSettings.t),CTsim.x2_sim)
    hold off
    grid on
    legend('x_1','x_2')
    title('Trajectories of x_1 and x_2 in CT simulation')
    xlabel('Time [s]')
    ylabel('Amplitude')
end

end