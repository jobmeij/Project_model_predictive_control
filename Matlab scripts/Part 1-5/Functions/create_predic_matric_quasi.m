%% MPC course
% create prediction matrices required for MPC
% create PHI and GAMMA
% inputs are A and B system matrices and the horizon N

function [QuasiLPV] = create_predic_matric_quasi(QuasiLPV)

A = QuasiLPV.Ad;
B = QuasiLPV.Bd;
N = QuasiLPV.N;

p = sym('p', [1 N]);
% x_in = sym(x_in)
n_states = 2;
n_inputs = 1;
% PHI = zeros(n_states*N,n_states);
% GAMMA = zeros(n_states*N,n_inputs*N);
% A_element = sym(2,2)
clear PHI
for i = 1:N
    if i == 1
        A_element = A(p(i));
    else
        A_element = A_element*A(p(i));
    end
    startrow = (i-1)*n_states+1;
    endrow = startrow+n_states-1;
    PHI(startrow:endrow, 1:n_states) = A_element;
%     A_element = A_element*A(x(i+1));
end    
% PHI_func = matlabFunction(PHI, 'vars', {x});
% x_in = [0 0 0 0 0];
% PHI_func(x_in)

element = sym('B', [2,1]);
element = subs(element,element,B);
clear GAMMA temp
for i = 1:N
    temp = sym('temp', [n_states*N, n_inputs*N]);
    temp = subs(temp ,temp ,zeros(n_states*N, n_inputs*N));
    amount = N-i+1;                            
    Ar = repmat(element, 1, amount);                                  
    Ac = mat2cell(Ar, n_states, repmat(n_inputs,1,amount));    
    Out = blkdiag(Ac{:});
    temp(((i-1)*n_states+1):end,1:amount*n_inputs) = Out;
    if i == 1
    	GAMMA = temp;
    else
        GAMMA = GAMMA+temp;
    end
    element = A(p(i))*element; 
end    

% GAMMA_func = matlabFunction(GAMMA, 'vars', {x});
% x_in = [1 0 0 0 0]
% GAMMA_func(x_in)
% GAMMA

% QuasiLPV.PHI = PHI_func;
% QuasiLPV.GAMMA = GAMMA_func;

QuasiLPV.Phi = PHI;
QuasiLPV.Gamma = GAMMA;