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

% AGENTS DEFINITION
NN = 5;

% min{z1...zN} sum{1,N}(ci'*zi)
% subj: 	sum{1,N} (Hi*zi)=b == sum{1,N}(Hi*zi-bi)=0
%       	zi € Pi, i€{1...N}

% Hi*zi-bi==gi(zi)
b=progenres.b; c=progenres.c; d=progenres.d; D=progenres.D;
H=progenres.H; LB=progenres.LB; UB=progenres.UB; %%Do we need LB,UB?
%we still need to implement the Pi containing the limits on Dizi<=di



%%
% CENTRALIZED SOLUTION
options = optimoptions('linprog','Display','none');
[~, fopt, exit_flag] = linprog(c,[],[],H,b,LB,UB,options); %%linprog(c,D,d,H,b,[],[],options); this should be the correct linprog

if exit_flag ~= 1
  fprintf(2,'A problem occurred in the centralized solution\n');
  return;
end

fprintf('Centralized optimal cost is %.4g\n',fopt);


%%
% WEIGHTED ADJACENCY MATRIX
% for averaging mechanism
p = 0.1;
[AA_NW, AA] = binomialGraph(1, NN, 'doubly');

%%

MAXITERS = 5e4;

ZZ = zeros(NN,MAXITERS);
ZZ_RA = zeros(NN,MAXITERS);%running average

LA = zeros(NN,MAXITERS);
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
  
  alpha_t = 0.1*(1/tt)^0.6; % 1/tt^alpha with alpha in (0.5, 1]

  for ii=1:NN
    N_ii = find(AA_NW(:,ii) == 1)';

    % update: averaging
    % vi =  sum_{k:1:N} aik*muk^t
    VV(ii) = AA(ii,ii)*LA(ii,tt);
    
    for jj = N_ii
      VV(ii) = VV(ii) + AA(ii,jj)*LA(jj,tt); %%update on v
    end
  end

% Primal Update: local minimization  
% zi^t+1=argmin{zi€Pi}(f(zi) + (vi^t+1)' * gi(zi))
%       =argmin{zi€Pi}( (ci'+ (vi^t+1)'*Hi)*zi)

%   for ii=1:NN
%      XX(ii,tt) = linprog(c(ii)+VV(ii)*AA_LP(ii),[],[],[],[],LB(ii),UB(ii),options);
%   end
%%DO WE NEED THIS? maybe just implement the check with Dz<=d
  for ii=1:NN
    if c(ii)+VV(ii)*H(ii)>=0
        ZZ(ii,tt) = LB(ii);
    else
        ZZ(ii,tt) = UB(ii);
    end
  end
  
  % Running average
  for ii=1:NN
      if tt==1
          ZZ_RA(ii,tt) = ZZ(ii,tt);
      else
          ZZ_RA(ii,tt) = (1/tt)*((tt-1)*ZZ_RA(ii,tt-1)+ZZ(ii,tt));
      end
  end
  
  % Dual Update
  
  % mui^t+1=max{0, vi^t+1 + alfa^t *(Hi*zi^t+1 -bi) }
  % mu^{t+1} = mu^t + gamma^t* ( sum_i grad_q_i(mu^t) )

  for ii=1:NN
    grad_ii = ZZ(ii,tt)-b(ii);

    LA(ii,tt+1) = VV(ii) + alpha_t*grad_ii;
  end
  
  % Performance check
	LA_avg = mean(LA(:,tt));

  for ii=1:NN
    ff_ii = c(ii)*ZZ(ii,tt);
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = c(ii)*ZZ_RA(ii,tt);
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = c(ii)*ZZ(ii,tt) + LA(ii,tt)*(ZZ(ii,tt)-b(ii));
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = c(ii)*ZZ_RA(ii,tt) + LA(ii,tt)*(ZZ_RA(ii,tt)-b(ii));
    dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;

    consensus_err(tt) = consensus_err(tt) + norm(LA(ii,tt) - LA_avg);
  end
end
%%
% This is computing only the last value of MAXITERS for all the agents [for plot]
tt = MAXITERS;
fprintf('Iteration n. %d\n',tt);
LA_avg = mean(LA(:,tt));

for ii=1:NN
  N_ii = find(AA_NW(:,ii) == 1)';
  
  % v_i =  sum_{j \in N_i U { i } } w_ij mu_j^t
  VV(ii) = AA(ii,ii)*LA(ii,tt);
  
  for jj = N_ii
    VV(ii) = VV(ii) + AA(ii,jj)*LA(jj,tt);
  end
end

for ii=1:NN
    if c(ii)+VV(ii)*H(ii)>=0
        ZZ(ii,tt) = LB(ii);
    else
        ZZ(ii,tt) = UB(ii);
    end
    ZZ_RA(ii,tt) = (1/tt)*((tt-1)*ZZ_RA(ii,tt-1)+ZZ(ii,tt));

    ff_ii = c(ii)*ZZ(ii,tt);
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = c(ii)*ZZ_RA(ii,tt);
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = c(ii)*ZZ(ii,tt) + LA(ii,tt)*(ZZ(ii,tt)-b(ii));
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = c(ii)*ZZ_RA(ii,tt) + LA(ii,tt)*(ZZ_RA(ii,tt)-b(ii));
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
  plot(1:MAXITERS,LA, 'LineWidth',2);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('\LA_i^t')

%%
figure
  plot(1:MAXITERS,sum(ZZ,1)-b, 'LineWidth',2);
  hold on, grid on, zoom on
  plot(1:MAXITERS,sum(ZZ_RA,1)-b, 'LineWidth',2);
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
