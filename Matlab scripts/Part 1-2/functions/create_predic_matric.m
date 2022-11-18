%% MPC course
% create prediction matrices required for MPC
% create PHI and GAMMA
% inputs are A and B system matrices and the horizon N

function [PHI, GAMMA] = create_predic_matric(LinMPC)

A = LinMPC.Systemd.A;
B = LinMPC.Systemd.B;
N = LinMPC.N;

n_states = size(A,1);
n_inputs = size(B,2);
PHI = zeros(n_states*N,n_states);
GAMMA = zeros(n_states*N,n_inputs*N);
for i = 1:N
    startrow = (i-1)*n_states+1;
    endrow = startrow+n_states-1;
    PHI(startrow:endrow, 1:n_states) = A^i;  
end    
element = B;
for i = 1:N
    temp = zeros(n_states*N,n_inputs*N);
    amount = N-i+1;                            
    Ar = repmat(element, 1, amount);                                  
    Ac = mat2cell(Ar, n_states, repmat(n_inputs,1,amount));    
    Out = blkdiag(Ac{:});
    temp(((i-1)*n_states+1):end,1:amount*n_inputs) = Out;
    GAMMA = GAMMA+temp;
    element = A*element; 
end    
% 
% PHI
% GAMMA
