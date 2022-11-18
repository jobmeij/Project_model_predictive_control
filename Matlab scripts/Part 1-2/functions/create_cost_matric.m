%% MPC course
% create prediction matrices required for MPC
% create Omega and Psi
% inputs are Q, R, P and the horizon N

function [Omega, Psi] = create_cost_matric(Q, R, P ,N)

Q, R, P ,N

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

Omega, Psi
LinMPC.Omega, LinMPC.Psi
