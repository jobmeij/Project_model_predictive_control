function [refMPC] = Construct_MPC_constraints_reftracking(refMPC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to calculate MPC constraints matrices L and W using E,M, D and c
% (MPC Course (5LMB0, lecture 1, slide 65)
% 
% Inputs are the constraint in the following form:
%   Input constraints: u_low <= u_i <= u_high
%   State constraints: x_low <= x_i <= x_high
%
% If the input(s) or state(s) are unconstraint leave inputs empty ([])
% if some of the input(s) or state(s) are unconstraint but not all of them
% use (-/+ inf) for the unconstraint state(s) of input(s) 
%
% The inputs Mn and bn are the terminal state constraints for the 
% invariant terminal set. If not desired, use []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Phi = refMPC.Phi;
Gamma = refMPC.Gamma;
u_low = refMPC.u_low;
u_high = refMPC.u_high;
x_low = refMPC.x_low;
x_high = refMPC.x_high;
N = refMPC.N;
Mn = [];
bn = [];


n_states = size(Phi, 2);
n_inputs = size(Gamma, 2)/N;

n_inputconstr_low = size(u_low, 1);
n_stateconstr_low = size(x_low, 1);
n_inputconstr_high = size(u_high, 1);
n_stateconstr_high = size(x_high, 1);

%TODO implement Mn and bn
n_terminal = size(Mn, 1);

if n_terminal == 0
    given_term_const = false;
else
    given_term_const = true;
end

% Check if the dimentions of the inputs are correct
if n_inputconstr_low == 0
    constrained_input = false;
else
    constrained_input = true;
    if ((n_inputconstr_low ~= n_inputconstr_high) || (n_inputconstr_high ~= n_inputs))
        error('Given input constraints and/or Gamma have the wrong dimensions')
    end
end
if n_stateconstr_low == 0
    constrained_state = false;
else
    constrained_state = true;
    if ((n_stateconstr_low ~= n_stateconstr_high) || (n_stateconstr_high ~= n_states))
        error('Given state constraints and/or Phi have the wrong dimensions')
    end
end
n_constraints = (n_stateconstr_low+n_inputconstr_low)*2;

if constrained_state
    % Create matrix to put in diagonal
    Mi = zeros(n_inputconstr_low*2, n_stateconstr_low);
    Mi = [Mi; -eye(n_states)];
    Mi = [Mi; eye(n_states)];
%     Mi = [Mi; ];
    
    % Create M matrix
    M_top = zeros(size(Mi,1),N*size(Mi,2));
    M_content = repmat(Mi, 1, N-1);                                  
    M_arg = mat2cell(M_content, size(Mi,1), repmat(size(Mi,2),1,N-1));  
    M_main = [blkdiag(M_arg{:}), zeros((size(Mi,1)*(N-1)),size(Mi,2))];
    if given_term_const == false
        M_n = [zeros(n_states, n_states*(N-1)), -eye(n_states); ...
            zeros(n_states, n_states*(N-1)), eye(n_states)];
    else
        M_n = [zeros(n_terminal, n_states*(N-1)), Mn];
    end
    
    M = [M_top; M_main; M_n];
    
    % Create D matrix
    if given_term_const == false
        D = [Mi; zeros(size(Mi,1)*(N-1), n_states);zeros(n_states*2,n_states)];
    else
        D = [Mi; zeros(size(Mi,1)*(N-1), n_states);zeros(n_terminal,n_states)];
    end
   
else
    M = zeros(n_constraints*(N)+(2*n_states), n_states*N);
    D = zeros(n_constraints*(N)+n_states, n_states);
end    

if constrained_input
    % Create matrix to put in diagonal
%     Epsiloni = zeros(n_stateconstr_low*2, n_inputconstr_low);
    Epsiloni = -eye(n_inputs);
    Epsiloni = [Epsiloni; eye(n_inputs)];
    Epsiloni = [Epsiloni; zeros(n_stateconstr_low*2, n_inputconstr_low)];
    
    
    % Create M matrix
    if given_term_const == false
        E_bottom = zeros(n_states*2,N*size(Epsiloni,2));
    else
        E_bottom = zeros(n_terminal,N*size(Epsiloni,2));
    end
    E_content = repmat(Epsiloni, 1, N);                
    E_arg = mat2cell(E_content, size(Epsiloni,1), repmat(size(Epsiloni,2),1,N));
    E_main = blkdiag(E_arg{:});
    E = [E_main; E_bottom];
    
else
    E = zeros((n_constraints*(N)+n_states), n_inputs*N);
end  

% Construct constraint value matrix c
bi = [-u_low; u_high; -x_low;x_high;];
c = repmat(bi, (N), 1);
if given_term_const == false
    c = [c; -x_low;x_high];
else
    c =[c; bn];
end
% Calculate outputs
L = M*Gamma+E;
W = -D-M*Phi;

% M
% D
% E
% c

refMPC.W = W;
refMPC.L = L;
refMPC.c = c; 