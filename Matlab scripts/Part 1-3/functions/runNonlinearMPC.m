%% Part 1.3: Closed loop simulation with nonlinear MPC controller
function [NonLinMPC] = runNonlinearMPC(SimSettings,Model,MPC)
% Set workspace to this function
options = simset('SrcWorkspace','current');
Ts = SimSettings.Ts;    % Load sample time in this function workspace

yalmip('clear')

% Model data
NonLinMPC.x0 = [1 0]';
NonLinMPC.A = [1 0.1; -0.1 (1.2-0.2*NonLinMPC.x0(1)^2)];
NonLinMPC.B = [0; 0.1];               % Input is sample time Ts * u on velocity
NonLinMPC.C = [1 0];                  % Output is position
NonLinMPC.nx = 2;                     % Number of states
NonLinMPC.nu = 1;                     % Number of inputs

% MPC settings
NonLinMPC.Q = MPC.Q;
NonLinMPC.R = MPC.R;
NonLinMPC.N = MPC.N;
[~,NonLinMPC.P,~] = dlqr(NonLinMPC.A,NonLinMPC.B,NonLinMPC.Q,NonLinMPC.R);  % Compute K and P using discrete LQR

% Constraints (state and input)
NonLinMPC.umin = Model.uMin;
NonLinMPC.umax = Model.uMax;
NonLinMPC.x1min = Model.pMin;
NonLinMPC.x1max = Model.pMax;
NonLinMPC.x2min = Model.pDotMin;
NonLinMPC.x2max = Model.pDotMax;

% Creating state and input sdpvars
NonLinMPC.u = sdpvar(repmat(NonLinMPC.nu,1,NonLinMPC.N),repmat(1,1,NonLinMPC.N));
NonLinMPC.x = sdpvar(repmat(NonLinMPC.nx,1,NonLinMPC.N+1),repmat(1,1,NonLinMPC.N+1));

% Compute nonlinear MPC cost and constraints for prediction horizon
NonLinMPC.constraints = [];
NonLinMPC.cost = 0;
for k = 1:NonLinMPC.N
    NonLinMPC.cost = NonLinMPC.cost + NonLinMPC.x{k}'*NonLinMPC.Q*NonLinMPC.x{k} + NonLinMPC.u{k}'*NonLinMPC.R*NonLinMPC.u{k};
    NonLinMPC.constraints = [NonLinMPC.constraints, NonLinMPC.x{k+1} == NonlinearDiscreteModel(NonLinMPC.x{k},NonLinMPC.u{k})];
    NonLinMPC.constraints = [NonLinMPC.constraints, NonLinMPC.umin <= NonLinMPC.u{k}<= NonLinMPC.umax, NonLinMPC.x1min <=NonLinMPC.x{k+1}(1)<= NonLinMPC.x1max, NonLinMPC.x2min <=NonLinMPC.x{k+1}(2)<= NonLinMPC.x2max];
end
% Compute terminal cost and constraints
NonLinMPC.cost = NonLinMPC.cost + NonLinMPC.x{NonLinMPC.N}'*NonLinMPC.P*NonLinMPC.x{NonLinMPC.N};
NonLinMPC.constraints = [NonLinMPC.constraints, NonLinMPC.x1min <= NonLinMPC.x{k+1}(1) <= NonLinMPC.x1max, NonLinMPC.x2min <= NonLinMPC.x{k+1}(2) <= NonLinMPC.x2max];

% Inputs and outputs
NonLinMPC.parameters_in = NonLinMPC.x{1};
NonLinMPC.solutions_out = [NonLinMPC.u{:}];

% Compute control output using optimizer (nonlinear optimzation solver)
NonLinMPC.controller = optimizer(NonLinMPC.constraints,NonLinMPC.cost,[],NonLinMPC.parameters_in,NonLinMPC.solutions_out);
x = SimSettings.x0;

if (SimSettings.NonLinMPCplot)
    clf;
    hold on
end

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

% Run nonlinear MPC simulation
for i = 1:(NonLinMPC.Length-1)
    [solutions,diagnostics] = NonLinMPC.controller{x};
    if diagnostics == 1
        error('Infeasible!');
    end
    NonLinMPC.U = solutions;                                    % Input prediction horizon

    % Plotting if enabled
    if (SimSettings.NonLinMPCplot)
        cla;
        hold on
        stairs(i:i+length(NonLinMPC.U)-1,NonLinMPC.U,'r')
        stairs(i:i+length(NonLinMPC.X)-1,NonLinMPC.X(1,:),'b')
        stairs(i:i+length(NonLinMPC.X)-1,NonLinMPC.X(2,:),'g')
        hold off
        grid on
        title('Predicted input U and states X_1 X_2')
        xlabel('Sample [N]')
        ylabel('Amplitude [-]')
        legend('U','X_1','X_2','Location','Best')
    end
    
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
    
    % Run simulink model
    sim('CTmodel.slx',[],options);
    
    % Save results
    NonLinMPC.x1simulink(i+1) = x1simu(end);
    NonLinMPC.x2simulink(i+1) = x2simu(end);

    x = [x1simu(end) x2simu(end)]';
   
    if (SimSettings.NonLinMPCplot)
        pause(0.025)        % Add a pause for plotting trajectory
    end
end

if (SimSettings.NonLinMPCplot)
    % Plotting trajectories of input and states
    figure()
    hold on
    stairs(1:NonLinMPC.Length,NonLinMPC.Utrajectory,'r')
    stairs(1:NonLinMPC.Length,NonLinMPC.Xtrajectory(:,1),'b')
    stairs(1:NonLinMPC.Length,NonLinMPC.Xtrajectory(:,2),'g')
    hold off
    title('Trajectories of u, x_1 and x_2')
    grid on
    xlabel('Sample [N]')
    ylabel('Aamplitude [-]')
    legend('u','x_1','x_2')
    
    % Plot trajectories of states on x and y axis
    figure()
    stairs(NonLinMPC.Xtrajectory(:,1),NonLinMPC.Xtrajectory(:,2))
    grid on
    title('X_1 and X_2')
    xlabel('X_1')
    ylabel('X_2')
    
    % Simulink simulation plot
    figure()
    hold on
    stairs(1:NonLinMPC.Length,NonLinMPC.usimulink)
    stairs(1:NonLinMPC.Length,NonLinMPC.x1simulink)
    stairs(1:NonLinMPC.Length,NonLinMPC.x2simulink)
    hold off
    grid on
    xlabel('Sample [N]')
    ylabel('Amplitude [-]')
    title('Simulink simulation result')
    legend('u','x_1','x_2','Location','Best')
end

end

% Nonlinear discrete-time model for nonlinear MPC
function [xplus] = NonlinearDiscreteModel(x,u)
    xplus = [x(1) + 0.1*x(2); -0.1*x(1) + (1.2-0.2*x(1)^2)*x(2) + 0.1*u];
end