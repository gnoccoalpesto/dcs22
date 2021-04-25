clear,clc,close all
rng(1);

% DCS 2020
% PROJECT 3
% GROUP 22: CANELLO
%           CERRI
%           RONCATO

%task1

% adds to path the folder containing this file and library subfolders
addpath(fileparts(which('task2dcs22')));
addpath(fullfile(fileparts(which('task2dcs22')),'/both'));
addpath(fullfile(fileparts(which('task2dcs22')),'/task1'));


%% 
disp('generating problem ...')

%number of agents
AgN=3;
%dimension of agents
Agni=4;

probdata=progen(AgN,Agni);

stocasticity='doubly';
graphtype='binomial';
graphprob=1;% €(0,1]
[ANW,A]=myGraph(AgN,stocasticity,graphtype,graphprob);
%ANW unused (since it's easily identifiable from A)

 
disp('generated!')

b=probdata.b;c=probdata.c;
d=probdata.d;D=probdata.D;
H=probdata.H;UB=probdata.UB;
LB=probdata.LB;

%%
disp('solving with dual subgradient method...')

[primal_cost,dual_cost,primal_cost_RA,dual_cost_RA,consensus_error,...
    lambda,ZZ,ZRA,maxIters] =dualsub(AgN,Agni,A,ANW,b,c,d,D,H,LB,UB);
%d,D,LB unused

%%
disp('centralized solution');

options = optimoptions('linprog','Display','none');
[~, fopt, exit_flag] = linprog(c,D,d,H,b,LB,UB,options);

% % if exit_flag ~= 1
% %   fprintf(2,'A problem occurred in the centralized solution\n');
% %   return;
% % end

fprintf('Centralized optimal cost is %.4g\n',fopt);

%%
disp('printing...')
% figure
%   semilogy(1:maxIters,abs(primal_cost(1:maxIters)-fopt), 'LineWidth',2);
%   hold on, grid on, zoom on
%   semilogy(1:maxIters,abs(primal_cost_RA(1:maxIters)-fopt), 'LineWidth',2);
%   semilogy(1:maxIters,abs(dual_cost(1:maxIters)-fopt), 'LineWidth',5);
%   semilogy(1:maxIters,abs(dual_cost_RA(1:maxIters)-fopt), 'LineWidth',2);
%   xlabel('t')
%   ylabel('cost error')
%   legend('primal cost','primal cost with running avg',...
%                 'dual cost','dual cost with running avg')
% 
% % %
% % hold on
% % for kk=1:Agni
% %     plot(1:maxIters,lambda(:,kk,:), 'LineWidth',2);
% %     hold on, grid on, zoom on
% %     xlabel('t')
% %     ylabel('\lambda_i^t')
% % end
% % hold off
% % %
% zum=reshape(sum(ZZ),Agni,maxIters);
% zumra=reshape(sum(ZRA),Agni,maxIters);
% for kk=1:Agni
%     figure
%       plot(1:maxIters,zum(kk,:)-b(kk), 'LineWidth',2);
%       hold on, grid on, zoom on
%       plot(1:maxIters,zumra(kk,:)-b(kk), 'LineWidth',2);
%       xlabel('t')
%       ylabel('x_1^t + ... + x_N^t - b')
%       legend('z','z from running avg')
% %       ylim([(min(min(zum,[],'all')-0.2*max(zum,[],'all'),...
% %           min(zumra,[],'all')-0.2*max(zumra,[],'all'))),...
% %           1.2*max(max(zum,[],'all'),max(zumra,[],'all'))]) 
% end  
consMatr4print=reshape(consensus_error,Agni,maxIters);
for kk=11:Agni%works on single iteraction
    figure
    consRow4print=consMatr4print(kk,:);
      semilogy(1:maxIters,consRow4print(1:maxIters), 'LineWidth',2);
%       hold on, grid on, zoom on
%       xlabel('t')
%       ylabel('consensus error')
end