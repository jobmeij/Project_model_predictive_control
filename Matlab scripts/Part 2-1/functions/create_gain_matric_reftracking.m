% Create gain matrix for MPC
function [refMPC] = create_gain_matric_reftracking(refMPC)

Psi = refMPC.Psi;
Gamma = refMPC.Gamma;
Omega = refMPC.Omega;
Phi = refMPC.Phi;
N = refMPC.N;


G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;

refMPC.G = G;
refMPC.F = F;
end
