% Create gain matrix for MPC


function [QuasiLPV] = create_gain_matric_quasi(QuasiLPV)

Psi = QuasiLPV.Psi;
Gamma = QuasiLPV.Gamma;
Omega = QuasiLPV.Omega;
Phi = QuasiLPV.Phi;
N = QuasiLPV.N;

G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;

% n_inputs = size(Gamma,2)/N;
% extrax_matrix = [eye(n_inputs), repmat(zeros(n_inputs,n_inputs),1,(N-1))];
% inv(G)*F;
% Kmpc = -extrax_matrix*inv(G)*F;Kmpc

% QuasiLPV.Kmpc = Kmpc;
QuasiLPV.G = G;
QuasiLPV.F = F;