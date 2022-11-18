% Create gain matrix for MPC


function [Kmpc, G, F] = create_gain_matric_reftracking(Psi, Gamma, Omega, Phi, N)

G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;

n_inputs = size(Gamma,2)/N;
extrax_matrix = [eye(n_inputs), repmat(zeros(n_inputs,n_inputs),1,(N-1))];
inv(G)*F;
Kmpc = -extrax_matrix*inv(G)*F;
