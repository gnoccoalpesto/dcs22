clear,clc,close all
rng(1);
%% 
disp('generating problem ...')

agN=3;
agM=4;

probdata=progen(agN,agM);

stocasticity='doubly';
graphtype='binomial';
graphprob=1;% €(0,1]
[ANW,A]=myGraph(agN,stocasticity,graphtype,graphprob);
%ANW unused (since it's easily identifiable from A)

 
disp('generated!')

b=probdata.b;c=probdata.c;
%d=probdata.d;D=probdata.D;
H=probdata.H;UB=probdata.UB;
LB=probdata.LB;

%%
disp('solving with dual subgradient method...')

[primal_cost,dual_cost,primal_cost_RA,dual_cost_RA,consensus_error,...
    lamba,ZZ,ZRA,maxIters] =dualsub(agN,A,b,c,[],[],H,[],UB);
%d,D,LB unused

%%
disp('centralized solution');

options = optimoptions('linprog','Display','none');
[~, fopt, exit_flag] = linprog(c,[],[],H,b,LB,UB,options);

% % if exit_flag ~= 1
% %   fprintf(2,'A problem occurred in the centralized solution\n');
% %   return;
% % end

fprintf('Centralized optimal cost is %.4g\n',fopt);

%% 
disp('printing...')

figure
  semilogy(1:maxIters,abs(primal_cost(1:maxIters)-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogy(1:maxIters,abs(primal_cost_RA(1:maxIters)-fopt), 'LineWidth',2);
  semilogy(1:maxIters,abs(dual_cost(1:maxIters)-fopt), 'LineWidth',2);
  semilogy(1:maxIters,abs(dual_cost_RA(1:maxIters)-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('cost error')
  legend('primal cost','primal cost with running avg','dual cost','dual cost with running avg')
  
% %
% figure
%   plot(1:MAXITERS,lam, 'LineWidth',2);
%   hold on, grid on, zoom on
%   xlabel('t')
%   ylabel('\mu_i^t')
% 
% %
% figure
%   plot(1:MAXITERS,sum(zz,1)-bb_centr, 'LineWidth',2);
%   hold on, grid on, zoom on
%   plot(1:MAXITERS,sum(XX_RA,1)-bb_centr, 'LineWidth',2);
%   xlabel('t')
%   ylabel('x_1^t + ... + x_N^t - b')
%   legend('x','x from running avg')
%   
%   
% %
% figure
%   semilogy(1:MAXITERS,consensus_err(1:MAXITERS), 'LineWidth',2);
%   hold on, grid on, zoom on
%   xlabel('t')
%   ylabel('consensus error')