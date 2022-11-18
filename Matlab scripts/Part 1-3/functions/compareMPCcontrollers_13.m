%% Part 1.3: Compare the nonlinear and linear MPC controllers
function [] = compareMPCcontrollers_13(~, LinMPC, NonLinMPC)

% Simulink simulation plot
figure()
hold on
plot(0:(NonLinMPC.Length-1),NonLinMPC.usimulink,'r')
plot(0:(NonLinMPC.Length-1),NonLinMPC.x1simulink,'g')
plot(0:(NonLinMPC.Length-1),NonLinMPC.x2simulink,'b')
plot(0:(LinMPC.Length-1),LinMPC.u_sim,'--')
plot(0:(LinMPC.Length-1),LinMPC.x1_sim,'--')
plot(0:(LinMPC.Length-1),LinMPC.x2_sim,'--')
hold off
grid on
xlabel('Sample [k]')
ylabel('Amplitude [-]')
title('Comparing linear and nonlinear MPC controllers')
legend('Nonlinear u','Nonlinear x_1','Nonlinear x_2','Linear u','Linear x_1','Linear x_2','Location','Best')

end