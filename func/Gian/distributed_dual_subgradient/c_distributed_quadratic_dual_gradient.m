%
% Matlab exercising
% Dual Gradient
% Prof. Giuseppe Notarstefano
%
% Code
% Dr. Ivano Notarnicola, Andrea Testa
% Bologna, 25/11/2020
%
clear all; clc; close all;
rng(2)

NN = 5;

% f_i(x_i) = 0.5*(x_i - a_i)^2
% grad_f_i(x_i) =  x_i - a_i

aa = rand(NN,1);
% bb = 5*rand;
bb = -5*rand;

%% Centralized Solution

% min_{x_1,...,x_N} f(x) = \sum f_i(x_i)
options_qp = optimoptions('quadprog','Display','none');
[xopt, fopt, exit_flag] = ...
         quadprog(eye(NN),-aa,ones(NN,1)',bb,[],[],[],[],[],options_qp);

for ii=1:NN
  fopt = fopt + 0.5*aa(ii)^2;
end

if exit_flag ~= 1
  fprintf(2,'A problem occurred in the centralized solution\n');
  return;
end

fprintf('Centralized optimal cost is %.4g\n',fopt);

%%
p = 0.1;
[AA_NW, AA] = binomialGraph(1, NN, 'doubly');
%%

MAXITERS = 3e3;

XX = zeros(NN,MAXITERS);

MU = zeros(NN,MAXITERS);
VV = zeros(NN,1);

primal_cost = zeros(MAXITERS,1);
dual_cost   = zeros(MAXITERS,1);

consensus_err = zeros(MAXITERS,1);

for tt = 1:MAXITERS-1
  if mod(tt,100)==0
      fprintf('Iteration n. %d\n',tt);
  end

  gamma_t = 0.1*(1/tt)^0.6; % 1/tt^alpha with alpha in (0.5, 1]

  for ii=1:NN
    N_ii = find(AA_NW(:,ii) == 1)';

    % v_i =  sum_{j \in N_i U { i } } w_ij mu_j^t
    VV(ii) = AA(ii,ii)*MU(ii,tt);
    
    for jj = N_ii
      VV(ii) = VV(ii) + AA(ii,jj)*MU(jj,tt);
    end
  end

  % Primal Update
  for ii=1:NN
    XX(ii,tt) = quadprog(1,-aa(ii)+VV(ii),[],[],[],[],[],[],[],options_qp);
  end
  
	% Dual Update
  % mu^{t+1} = mu^t + gamma^t* ( sum_i grad_q_i(mu^t) )
  for ii=1:NN
    grad_ii = XX(ii,tt)-bb/NN;

    MU(ii,tt+1) = VV(ii) + gamma_t*grad_ii;
    MU(ii,tt+1) = max(MU(ii,tt+1),0);    
  end
  
  % Performance check
	MU_avg = mean(MU(:,tt));

  for ii=1:NN
    ff_ii = 0.5*(XX(ii,tt) - aa(ii))^2;
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    qq_ii = ff_ii + MU(tt)*(XX(ii,tt)-bb/NN);
    dual_cost(tt) = dual_cost(tt) + qq_ii;

    consensus_err(tt) = consensus_err(tt) + norm(MU(ii,tt) - MU_avg);
  end
end
%%
% Last value [for plot]
tt = MAXITERS;
MU_avg = mean(MU(:,tt));

for ii=1:NN
  N_ii = find(AA_NW(:,ii) == 1)';
  
  % v_i =  sum_{j \in N_i U { i } } w_ij mu_j^t
  VV(ii) = AA(ii,ii)*MU(ii,tt);
  
  for jj = N_ii
    VV(ii) = VV(ii) + AA(ii,jj)*MU(jj,tt);
  end
end

for ii=1:NN
    XX(ii,tt) = quadprog(1,-aa(ii)+VV(ii),[],[],[],[],[],[],[],options_qp);

    ff_ii = 0.5*(XX(ii,tt) - aa(ii))^2;
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    qq_ii = ff_ii + MU(tt)*(XX(ii,tt)-bb/NN);
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    consensus_err(tt) = consensus_err(tt) + norm(MU(ii,tt) - MU_avg);
end

%%
figure
  semilogy(1:MAXITERS,abs(primal_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogy(1:MAXITERS,abs(dual_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('cost error')
  legend('primal cost','dual cost')
  
%%
figure
  plot(1:MAXITERS,MU, 'LineWidth',2);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('\mu_i^t')

%%
figure
  plot(1:MAXITERS,sum(XX,1)-bb, 'LineWidth',2);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('x_1^t + ... + x_N^t - b')
  
%%
figure
  semilogy(1:MAXITERS,consensus_err(1:MAXITERS), 'LineWidth',2);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('consensus error')