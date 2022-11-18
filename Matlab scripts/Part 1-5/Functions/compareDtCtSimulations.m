%% Part 1.3: Comparing discrete- and continuous time simulations
function [] = compareDtCtSimulations(DTsim,CTsim,SimSettings)
% 1D plotting
figure()
set(gcf,'Position',[10 100 1200 400])
%
subplot 121 % CT system
hold on
plot(CTsim.t_sim,CTsim.u_sim)
plot(CTsim.t_sim,CTsim.x1_sim)
plot(CTsim.t_sim,CTsim.x2_sim)
hold off
grid on
xlabel('Time [s]', 'Interpreter','latex')
ylabel('Amplitude [-]', 'Interpreter','latex')
legend('$u_{mpc}$','$\dot{p}(t)$','$\ddot{p}(t)$','Location','Best','Interpreter','latex')
title('Continuous-time simulation', 'Interpreter','latex')
%
subplot 122 % DT system
hold on
stairs(0:length(SimSettings.t)-1,DTsim.u)
stairs(0:length(SimSettings.t)-1,DTsim.x1)
stairs(0:length(SimSettings.t)-1,DTsim.x2)
hold off
grid on
xlabel('Sample [k]', 'Interpreter','latex')
ylabel('Amplitude [-]', 'Interpreter','latex')
legend('$u_{mpc}$','$x_{1}(k)$','$x_{2}(k)$','Location','Best','Interpreter','latex')
title('Discrete-time simulation', 'Interpreter','latex')

% 2D plotting
figure()
set(gcf,'Position',[100 100 600 400])
%
hold on
plot(CTsim.x1_sim,CTsim.x2_sim)
plot(DTsim.x1,DTsim.x2)
plot(SimSettings.x0(1),SimSettings.x0(2),'.','Color','Black','MarkerSize',20)
hold off
grid on
title('Continuous-time vs. discrete-time simulation', 'Interpreter','latex')
xlabel('$\dot{p}(t), x_{1}(k)$', 'Interpreter','latex')
ylabel('$\ddot{p}(t), x_{2}(k)$', 'Interpreter','latex')
legend('Continuous-time','Discrete-time','Initial condition','Location','Best','Interpreter','latex')

end
