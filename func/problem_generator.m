clear all; close all; clc;

N = 3; % NUMBER OF AGENTS 
M = 4; % SIZE OF AGENTS' LOCAL OPTIM. VAR. ATTENTION! IT MUST STAND M>N

%%% MODEL A
% % profit associated with assigning job j to machine i
P_ = randi([10 25],N,M);
tt = max(P_,[],'all')+1;
P  = tt - P_;

% claim on the capacity of machine i by job j
W = randi([5 25],N,M);
C = zeros(N,1);

for ijk = 1:N 
    ijk
    J_i = [];
    for jjj = 1:M
      p_j = P(:,jjj);
      [~, ind_m] = min(p_j);
      ind_ = find(ind_m == ijk);
      if ~isempty(ind_)
        J_i = [J_i jjj];
      end
    end
    c_v = zeros(N,1);
    if ~isempty(J_i)
        for lll = 1:N
            w_l = W(lll,:);
            for jkl=1:numel(J_i)
                c_v(lll)=c_v(lll)+w_l(J_i(jkl));
            end
        end
    end
    C(ijk) = 9*(M/N) + 0.4*max(c_v); 
end



A_eq = repmat(eye(M,M), [1,N]); % == [...Hj...]
b_eq = ones(M,1); % == b

% machine capacity, local constraints
A = zeros(N,M*N);
for i = 1:N
    A(i,M*(i-1)+1:M*(i-1)+M) = W(i,:);
end
b = C;

%maximize profit
c = reshape(P',[1,M*N]);
LB = zeros((N*M),1);
UB = ones((N*M),1);

[xopt, fopt] = linprog(c,A,b,A_eq,b_eq,LB,UB)

% min{x} f' x | Ax=<b & Aeq x = beq & lb=<x=<ub

%clearvars -except c A b A_eq b_eq LB UB xopt fopt
