% Simulate continious time non-linear model for one timestep (Ts)
function [x_out] = sim_nonlinear_model(x_in, u_in, Ts_in)

global x umpc x1simu x2simu Ts

assignin('base','Ts',Ts_in)
assignin('base','x',x_in)
assignin('base','umpc',u_in)
% Ts = Ts_in;
% x = x_in;
% umpc = u_in;

% x_sim_in_1 = x_sim_in(1);
% x_sim_in_2 = x_sim_in(2);
% u_sim_in_u = u_sim_in;
% Ts_sim = Ts_sim_in;

% set_param('CT_NL_model/x1', 'Value', num2str(x_sim_in_1))
% set_param('CT_NL_model/x2', 'value', num2str(x_sim_in_2))
% set_param('CT_NL_model/u', 'value', num2str(u_sim_in_u))
% set_param('CT_NL_model', 'StopTime', num2str(Ts_sim))


sim('CTmodel.slx');

% u_sim(k) = umpc;
x_out(1) = x1simu(end);    % Select end value
x_out(2) = x2simu(end);    % Select end value