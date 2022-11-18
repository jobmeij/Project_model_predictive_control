% Simulate continious time non-linear model for one timestep (Ts)
function [x_out] = sim_nonlinear_model(x_in, u_in, Ts_in)

global x umpc x1simu x2simu Ts

assignin('base','Ts',Ts_in)
assignin('base','x',x_in)
assignin('base','umpc',u_in)

sim('CTmodel.slx');

% u_sim(k) = umpc;
x_out(1) = x1simu(end);    % Select end value
x_out(2) = x2simu(end);    % Select end value