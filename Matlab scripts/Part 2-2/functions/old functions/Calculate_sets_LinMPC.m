function [LinMPC] = Calculate_sets_LinMPC(LinMPC, figurenr)
% if figurenr is number -> variant set is plotted in that figure
% if figurenr = [] -> no figure is made
% if figurenr = 0 -> new figure is made

disp('Start calculating sets and linear cost matrices for linear MPC')

% Calculate terminal set and controller
[LinMPC.P, ~, K] = dare(LinMPC.Systemd.A, LinMPC.Systemd.B, LinMPC.Q, LinMPC.R);
LinMPC.K = -K;
LinMPC.Acl = LinMPC.Systemd.A+LinMPC.Systemd.B*LinMPC.K;


% % Calculate the invariant set
% disp('Calculate invariant set')
% LinMPC.LTImodel = LTISystem('A', LinMPC.Acl); %Acl
% LinMPC.Xset = Polyhedron([-eye(size(LinMPC.x_low,1)); eye(size(LinMPC.x_low,1))], [-LinMPC.x_low;LinMPC.x_high]);
% LinMPC.XUset = Polyhedron([-eye(size(LinMPC.u_low,1));eye(size(LinMPC.u_low,1))]*LinMPC.K,[-LinMPC.u_low; LinMPC.u_high]) & LinMPC.Xset;
% LinMPC.InvSetXU = LinMPC.LTImodel.invariantSet('X',LinMPC.XUset);

% % Calculate control admissable invariant set
% disp('Calculate control admissable set')
% model_CI = LTISystem('A', LinMPC.Systemd.A, 'B', LinMPC.Systemd.B);
% model_CI.x.min = LinMPC.x_low;
% model_CI.x.max = LinMPC.x_high;
% model_CI.u.min = LinMPC.u_low;
% model_CI.u.max = LinMPC.u_high;
% LinMPC.ContAdmInvSet = model_CI.invariantSet();

% Create MPC matrices
disp('Create cost matrices')
% LinMPC.Mn = LinMPC.InvSetXU.A;
% LinMPC.bn = LinMPC.InvSetXU.b;
[LinMPC.Phi, LinMPC.Gamma] = create_predic_matric(LinMPC.Systemd.A, LinMPC.Systemd.B, LinMPC.N);
[LinMPC.Omega, LinMPC.Psi] = create_cost_matric(LinMPC.Q, LinMPC.R, LinMPC.P ,LinMPC.N);
[LinMPC.Kmpc, LinMPC.G, LinMPC.F] = create_gain_matric(LinMPC.Psi, LinMPC.Gamma, LinMPC.Omega, LinMPC.Phi, LinMPC.N);
[LinMPC.W, LinMPC.L, LinMPC.c ] = Construct_MPC_constraints_with_terminal(LinMPC.Phi, LinMPC.Gamma, LinMPC.u_low, LinMPC.u_high, LinMPC.x_low, LinMPC.x_high, LinMPC.N, LinMPC.Mn, LinMPC.bn);

% % Calculate feasable set
% disp('Calculate feasable set')
% LinMPC.FeasSet = Polyhedron('A', [-LinMPC.W LinMPC.L], 'B', LinMPC.c);
% LinMPC.Feasiblesetx = projection(LinMPC.FeasSet,1:2);

% Plots
% if not(isempty(figurenr))
%     if figurenr == 0
%         figure() 
%     else 
%         figure(figurenr) 
%     end
%     % Plot control admissable invariant set
%     plot(LinMPC.ContAdmInvSet, 'Color','Blue', 'Alpha', 0.5);
%     hold on
%     % plot the feasable set
%     plot(LinMPC.Feasiblesetx, 'Color', 'Red', 'Alpha', 0.5); 
%     % plot the invariant set
%     plot(LinMPC.InvSetXU,'Color', 'yellow', 'Alpha', 0.5); 
%     title('State space sets', 'Interpreter','latex')
%     xlabel('State 1 (p)', 'Interpreter','latex')
%     ylabel('State 2 ($\dot{p}$)', 'Interpreter','latex')
%     legend('Control adm. inv. set','Feasable set (N=50)','Invariant set', 'location','Northwest')
%     hold off
% end
disp('Finished calculating sets and linear cost matrices for linear MPC')
