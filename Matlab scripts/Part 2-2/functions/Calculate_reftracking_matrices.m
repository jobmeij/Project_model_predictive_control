function [refMPC] = Calculate_reftracking_matrices(refMPC)

A = refMPC.Ad;
B = refMPC.Bd;
C = refMPC.Cd;
N = refMPC.N;
Gamma = refMPC.Gamma;
Omega = refMPC.Omega;
Psi = refMPC.Psi;


YA_content = repmat(A, 1, N);
YA_arg = mat2cell(YA_content, 2, repmat(2,1,N));
YA_main = [blkdiag(YA_arg{:}), zeros((2*(N)),2)];
YI_content = repmat(-eye(2), 1, N);
YI_arg = mat2cell(YI_content, 2, repmat(2,1,N));
YI_main = [zeros((2*(N)),2), blkdiag(YI_arg{:})];
YAI = [zeros(2,N*2) A-eye(2)];
YUL = [YA_main+YI_main; YAI];

YB_content = repmat(B, 1, N+1);
YB_arg = mat2cell(YB_content, 2, repmat(1,1,N+1));
YB_main = blkdiag(YB_arg{:});

YC_content = repmat(C, 1, N+1);
YC_arg = mat2cell(YC_content, 1, repmat(2,1,N+1));
YC_main = blkdiag(YC_arg{:});

YD = zeros(N+1, N+1);

SH_pre = [YUL YB_main; YC_main YD];
SH = inv(SH_pre);

H = SH(:,((2*(N+1))+1):end);

refMPC.Y = -2*[Gamma'*Omega Psi]*[zeros(3*N,2) eye(3*N) zeros(3*N,1)]*H;