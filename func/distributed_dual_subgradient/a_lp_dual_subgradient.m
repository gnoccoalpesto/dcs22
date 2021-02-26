clear all;close all;clc
rng(1)
NN = 5;

% min{z1...zN} sum{1,N}(ci'*zi)
% subj: 	sum{1,N} (Hi*zi)=b == sum{1,N}(Hi*zi-bi)=0
%       	zi € Pi, i€{1...N}

% Hi*zi-bi==gi(zi)



%weights
a = -1;   b = -10;
cc_LP = zeros(NN,1); % ci 
for ii = 1:NN
  cc_LP(ii) = (b-a)*rand(1 ,1) + a; % entries unif. in [a,b]
end

%declare Di and Gi here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
AA_LP = ones(1,NN);% Hi
bb_LP = zeros(NN,1);% b

%# chosing b s.t. solution always feasible
ll=randi([1,NN-1]);
a = 3*NN;  b = 4*NN;
for ii=1:NN
    bb_LP(ii) = (1/NN)*(LB(ll+1)+(b-a)*rand(1)+a -LB(ll));
end

bb_centr = sum(bb_LP);

%%
% Centralized Solution

% centralized (parallel) dual subgradient:
% zi^t+1= argmin{zi € Pi}(ci'*zi +(mu^t)'*(Hi*zi-bi)
%       = argmin{zi € Pi}( ci'+(mu^t)'*Hi)*zi -(mu^t)'*bi)
% mu^t+1 = max {0, mu^t+alfa^t *sum{1,N}(Hi*zi^t+1 -bi) }


options = optimoptions('linprog','Display','none');
[~, fopt, exit_flag] = linprog(cc_LP,AA_LP,bb_centr,[],[],LB,UB,options);

if exit_flag ~= 1
  fprintf(2,'A problem occurred in the centralized solution\n');
  return;
end

fprintf('Centralized optimal cost is %.4g\n',fopt);

%%

MAXITERS = 5e4;

XX = zeros(NN,MAXITERS);
XX_RA = zeros(NN,MAXITERS);

MU = zeros(MAXITERS,1);

primal_cost = zeros(MAXITERS,1);
dual_cost   = zeros(MAXITERS,1);
primal_cost_RA = zeros(MAXITERS,1);
dual_cost_RA   = zeros(MAXITERS,1);

tt_ra = 1;
NMIN=1;

for tt = 1:MAXITERS-1
  if mod(tt,100)==0
      fprintf('Iteration n. %d\n',tt);
  end
  % Primal Update
%   for ii=1:NN
%      XX(ii,tt) = linprog(cc_LP(ii)+MU(tt)*AA_LP(ii),[],[],[],[],LB(ii),UB(ii),options);
%   end
  for ii=1:NN
    if cc_LP(ii)+MU(tt)*AA_LP(ii)>=0
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
%   
  % Dual Update
  % mu^{t+1} = mu^t + gamma^t* ( sum_i grad_q_i(mu^t) )
   gamma_t = 0.1*(1/tt)^0.6; % 1/tt^alpha with alpha in (0.5, 1]

  
  MU(tt+1) = MU(tt);
  for ii=1:NN
    grad_ii = XX(ii,tt)-bb_LP(ii);
    
    MU(tt+1) = MU(tt+1) + gamma_t*grad_ii;
  end
	MU(tt+1) = max(0,MU(tt+1));
  
  % Performance check

  for ii=1:NN
    ff_ii = cc_LP(ii)*XX(ii,tt);
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = cc_LP(ii)*XX_RA(ii,tt);
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = cc_LP(ii)*XX(ii,tt) + MU(tt)*(XX(ii,tt)-bb_LP(ii));
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = cc_LP(ii)*XX_RA(ii,tt) + MU(tt)*(XX_RA(ii,tt)-bb_LP(ii));
    dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;
  end
end
%%
% Last value [for plot]
tt = MAXITERS;
fprintf('Iteration n. %d\n',tt);
for ii=1:NN   
    XX(ii,tt) = linprog(cc_LP(ii)+MU(tt)*AA_LP(ii),[],[],[],[],LB(ii),UB(ii),options);
    XX_RA(ii,tt) = (1/tt)*((tt-1)*XX_RA(ii,tt-1)+XX(ii,tt));
    ff_ii = cc_LP(ii)*XX(ii,tt);
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = cc_LP(ii)*XX_RA(ii,tt);
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = cc_LP(ii)*XX(ii,tt) + MU(tt)*(XX(ii,tt)-bb_LP(ii));
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = cc_LP(ii)*XX_RA(ii,tt) + MU(tt)*(XX_RA(ii,tt)-bb_LP(ii));
    dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;
end

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

% return
%% dual function 
% (remove the return)
mu = 0:0.01:10;
dual_=zeros(numel(mu),1);

for mm = 1:numel(mu)
        [~, fopt, ~] = linprog(cc_LP+mu(mm)*ones(NN,1),[],[],[],[],LB,UB,options);
        dual_(mm) = fopt-mu(mm)*bb_centr;
end
    
figure
plot(mu,dual_);grid on;zoom on;
xlabel('\mu')
ylabel('q(\mu)')



%% comments
% xi^t (blue one) oscillates, not converging
% running avarage converges 
%       XX_RA(ii,tt) = (1/tt)*((tt-1)*XX_RA(ii,tt-1)+XX(ii,tt))
% mu^t converging as expected for the dual update