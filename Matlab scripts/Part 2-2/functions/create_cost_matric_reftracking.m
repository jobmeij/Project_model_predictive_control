%% MPC course
% create prediction matrices required for MPC
% create Omega and Psi
% inputs are Q, R, P and the horizon N

function [refMPC] = create_cost_matric_reftracking(refMPC)

Q = refMPC.Q; 
R = refMPC.R; 
P = refMPC.P;
N = refMPC.N;

% Psi
[heigth, width] = size(R);
Ar = repmat(R, 1, N);
Ac = mat2cell(Ar, width, repmat(heigth,1,N));
Psi = blkdiag(Ac{:});

% Omega
[width, heigth] = size(Q);
Ar = repmat(Q, 1, N);
Ac = mat2cell(Ar, width, repmat(heigth,1,N));
Omega = blkdiag(Ac{:});
% Replace last matrix in omega with P
Omega((end-heigth+1):end,(end-width+1):end) = P;

refMPC.Omega = Omega;
refMPC.Psi = Psi;
