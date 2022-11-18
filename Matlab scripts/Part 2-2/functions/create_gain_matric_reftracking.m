% Create gain matrix for MPC


function [refMPC] = create_gain_matric_reftracking(refMPC)

Psi = refMPC.Psi;
Gamma = refMPC.Gamma;
Omega = refMPC.Omega;
Phi = refMPC.Phi;
N = refMPC.N;


G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;

% n_inputs = size(Gamma,2)/N;
% extrax_matrix = [eye(n_inputs), repmat(zeros(n_inputs,n_inputs),1,(N-1))];
% % Kmpc = -extrax_matrix*inv(G)*F;


refMPC.G = G;
refMPC.F = F;
