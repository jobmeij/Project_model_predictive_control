%% Part 1.3: Simulating the discretized model with given initial condition and zero input
function [DTsim] = runDiscreteTimeModel(SimSettings)
% Pre-allocate
DTsim.x1p = zeros(length(SimSettings.t),1);
DTsim.x2p = zeros(length(SimSettings.t),1);
DTsim.x1 = zeros(length(SimSettings.t),1);
DTsim.x1(1) = SimSettings.x0(1);
DTsim.x2 = zeros(length(SimSettings.t),1);
DTsim.x2(1) = SimSettings.x0(2);
DTsim.u = zeros(length(SimSettings.t),1);

% Simulate nonlinear discrete-time model in open loop:
%   x1(k+1) = x2(k)
%   x2(k+1) = (-0.1*x1(k)) + (1.2-(0.2*x1^2))*x2(k) + 0.1*u(k)
for k = 1:length(SimSettings.t)
    % Compute states
    DTsim.x1p(k) = DTsim.x1(k) + 0.1*DTsim.x2(k);
    DTsim.x2p(k) = -0.1*DTsim.x1(k) + (1.2-(0.2*DTsim.x1(k)^2))*DTsim.x2(k) + 0.1*DTsim.u(k);
    
    % Set x(k+1) to x(k) for next iteration 
    if (k ~= length(SimSettings.t))
        DTsim.x1(k+1) = DTsim.x1p(k);
        DTsim.x2(k+1) = DTsim.x2p(k);
    end
end

% Plotting of nonlinear DT model (if enabled)
if (SimSettings.DTplot)
    figure()
    hold on
    stairs(1:length(SimSettings.t),DTsim.x1)
    stairs(1:length(SimSettings.t),DTsim.x2)
    hold off
    grid on
    legend('x_1','x_2')
    title('Trajectories of x_1 and x_2 in DT simulation')
    xlabel('Sample [N]')
    ylabel('Amplitude')
end

end