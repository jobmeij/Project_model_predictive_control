%% Part 1.4: Compare the controllers
function [] = compareMPCcontrollers_14(LinMPC, NonLinMPC, QuasiLPV)

% Simulink simulation plot
figure()
hold on
plot(0:(NonLinMPC.Length-1),NonLinMPC.usimulink,'r')
plot(0:(NonLinMPC.Length-1),NonLinMPC.x1simulink,'g')
plot(0:(NonLinMPC.Length-1),NonLinMPC.x2simulink,'b')
plot(0:(LinMPC.Length-1),LinMPC.u_sim,'--r')
plot(0:(LinMPC.Length-1),LinMPC.x1_sim,'--g')
plot(0:(LinMPC.Length-1),LinMPC.x2_sim,'--b')
plot(0:(length(QuasiLPV.u_sol)-1),QuasiLPV.u_sol,':r')
plot(0:(length(QuasiLPV.u_sol)-1),QuasiLPV.x_sol(1,1:(end-1)),':g')
plot(0:(length(QuasiLPV.u_sol)-1),QuasiLPV.x_sol(2,1:(end-1)),':b')
hold off
grid on
xlabel('Sample [k]')
ylabel('Amplitude [-]')
title('Comparing quasiLPV, linear MPC and nonlinear MPC controllers')
legend('Nonlinear u','Nonlinear x_1','Nonlinear x_2','Linear u','Linear x_1','Linear x_2','QuasiLPV u','QuasiLPV x_1','QuasiLPV x_2','Location','NorthEast')

end