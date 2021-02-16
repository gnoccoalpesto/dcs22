%
% Matlab exercising
% Distributed Dual Gradient
% Prof. Giuseppe Notarstefano
%
% Code
% Dr. Ivano Notarnicola, Andrea Testa
% Bologna, 25/11/2020
%
clear all; clc; close all;
rng(1)
NN = 5;

% min{z1...zN} sum{1,N}(ci'*zi)
% subj: 	sum{1,N} (Hi*zi)=b == sum{1,N}(Hi*zi-bi)=0
%       	zi € Pi, i€{1...N}

% Hi*zi-bi==gi(zi)

% weights ci
a = -1;   b = -10;
cc_LP = zeros(NN,1);
for ii = 1:NN
  cc_LP(ii) = (b-a)*rand(1 ,1) + a; % entries unif. in [a,b]
end


% Pi definition
% Local Box Constraint X_i
LB = zeros(NN,1);
UB = zeros(NN,1);

a = 0;  b = 5;
c =   6;  d = 15;
for ii=1:NN
  LB(ii) = (b-a)*rand(1) + a; % lower bound in [a, b]
  UB(ii) = (d-c)*rand(1) + c; % upper bound in [c, d]
end

% Linear Coupling Constraint
AA_LP = ones(1,NN); % Hi
bb_LP = zeros(NN,1); % b

ll=randi([1,NN-1]);
a = 3*NN;  b = 4*NN;
for ii=1:NN
    bb_LP(ii) = (1/NN)*(LB(ll+1)+(b-a)*rand(1)+a -LB(ll));
end

bb_centr = sum(bb_LP);

%%
% Centralized Solution
options = optimoptions('linprog','Display','none');
[~, fopt, exit_flag] = linprog(cc_LP,AA_LP,bb_centr,[],[],LB,UB,options);

if exit_flag ~= 1
  fprintf(2,'A problem occurred in the centralized solution\n');
  return;
end

fprintf('Centralized optimal cost is %.4g\n',fopt);


%%
% weighted adjacency matrix for averaging mechanism
p = 0.1;
[AA_NW, AA] = binomialGraph(1, NN, 'doubly');
%%

MAXITERS = 5e4;

XX = zeros(NN,MAXITERS);
XX_RA = zeros(NN,MAXITERS);

MU = zeros(NN,MAXITERS);
VV = zeros(NN,1);

primal_cost = zeros(MAXITERS,1);
dual_cost   = zeros(MAXITERS,1);
primal_cost_RA = zeros(MAXITERS,1);
dual_cost_RA   = zeros(MAXITERS,1);

consensus_err = zeros(MAXITERS,1);

for tt = 1:MAXITERS-1
  if mod(tt,100)==0
      fprintf('Iteration n. %d\n',tt);
  end
  
  gamma_t = 0.1*(1/tt)^0.6; % 1/tt^alpha with alpha in (0.5, 1]

  for ii=1:NN
    N_ii = find(AA_NW(:,ii) == 1)';

    % update: averaging
    % vi =  sum_{k:1:N} aik*muk^t
    VV(ii) = AA(ii,ii)*MU(ii,tt);
    
    for jj = N_ii
      VV(ii) = VV(ii) + AA(ii,jj)*MU(jj,tt);
    end
  end

% Primal Update: local minimization  
% zi^t+1=argmin{zi€Pi}(f(zi) + (vi^t+1)' * gi(zi))
%       =argmin{zi€Pi}( (ci'+ (vi^t+1)'*Hi)*zi)

%   for ii=1:NN
%      XX(ii,tt) = linprog(cc_LP(ii)+VV(ii)*AA_LP(ii),[],[],[],[],LB(ii),UB(ii),options);
%   end
  for ii=1:NN
    if cc_LP(ii)+VV(ii)*AA_LP(ii)>=0
        XX(ii,tt) = LB(ii);
    else
        XX(ii,tt) = UB(ii);
    end
  end
  
  % Running average
  for ii=1:NN
      if tt==1
          XX_RA(ii,tt) = XX(ii,tt);
      else
          XX_RA(ii,tt) = (1/tt)*((tt-1)*XX_RA(ii,tt-1)+XX(ii,tt));
      end
  end
  
  % Dual Update
  
  % mui^t+1=max{0, vi^t+1 + alfa^t *(Hi*zi^t+1 -bi) }
  % mu^{t+1} = mu^t + gamma^t* ( sum_i grad_q_i(mu^t) )
  for ii=1:NN
    grad_ii = XX(ii,tt)-bb_LP(ii);

    MU(ii,tt+1) = VV(ii) + gamma_t*grad_ii;
    MU(ii,tt+1) = max(MU(ii,tt+1),0);    
  end
  
  % Performance check
	MU_avg = mean(MU(:,tt));

  for ii=1:NN
    ff_ii = cc_LP(ii)*XX(ii,tt);
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = cc_LP(ii)*XX_RA(ii,tt);
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = cc_LP(ii)*XX(ii,tt) + MU(ii,tt)*(XX(ii,tt)-bb_LP(ii));
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = cc_LP(ii)*XX_RA(ii,tt) + MU(ii,tt)*(XX_RA(ii,tt)-bb_LP(ii));
    dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;

    consensus_err(tt) = consensus_err(tt) + norm(MU(ii,tt) - MU_avg);
  end
end
%%
% Last value [for plot]
tt = MAXITERS;
fprintf('Iteration n. %d\n',tt);
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
    if cc_LP(ii)+VV(ii)*AA_LP(ii)>=0
        XX(ii,tt) = LB(ii);
    else
        XX(ii,tt) = UB(ii);
    end
    XX_RA(ii,tt) = (1/tt)*((tt-1)*XX_RA(ii,tt-1)+XX(ii,tt));

    ff_ii = cc_LP(ii)*XX(ii,tt);
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = cc_LP(ii)*XX_RA(ii,tt);
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = cc_LP(ii)*XX(ii,tt) + MU(ii,tt)*(XX(ii,tt)-bb_LP(ii));
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = cc_LP(ii)*XX_RA(ii,tt) + MU(ii,tt)*(XX_RA(ii,tt)-bb_LP(ii));
    dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;
end

%%
figure
  semilogy(1:MAXITERS,abs(primal_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogy(1:MAXITERS,abs(primal_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
  semilogy(1:MAXITERS,abs(dual_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  semilogy(1:MAXITERS,abs(dual_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('cost error')
  legend('primal cost','primal cost with running avg','dual cost','dual cost with running avg')
  
%%
figure
  plot(1:MAXITERS,MU, 'LineWidth',2);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('\mu_i^t')

%%
figure
  plot(1:MAXITERS,sum(XX,1)-bb_centr, 'LineWidth',2);
  hold on, grid on, zoom on
  plot(1:MAXITERS,sum(XX_RA,1)-bb_centr, 'LineWidth',2);
  xlabel('t')
  ylabel('x_1^t + ... + x_N^t - b')
  legend('x','x from running avg')
  
  
%%
figure
  semilogy(1:MAXITERS,consensus_err(1:MAXITERS), 'LineWidth',2);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('consensus error')
  
  %% conclusions
  % dual cost: converging rapidly (orange), then oscillate
  % (violet) update on running average
  % (blue plot) zi^t non converging, oscillating
  % (red) running average is converging
  %
  % convergence on the mu variables
  % convergence of running averages
