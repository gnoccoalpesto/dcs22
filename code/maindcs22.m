% DCS 2020
% PROJECT 3
% GROUP 22: CANELLO
%           CERRI
%           RONCATO

%% TASK 1
% 
%  this code is intended to solve a constaint coupled linea problem using the 
% distibuted dual subgradent algorithm.
%  Data is generated randomly, while its dimension is decided: a script
% generates the system graph, while another generates problem data.
%     

clear,clc,close all
rng(1);

% adds to path the folder containing this file and library subfolders
addpath(fileparts(which('task2dcs22')));
%NOTE: scriptname must coincide, mfilename can prove usefull in case this
% is executed from matlab command window
% adds function folders
addpath(fullfile(fileparts(which('task2dcs22')),'/both'));% general
addpath(fullfile(fileparts(which('task2dcs22')),'/task1'));% specific for T1


%% data generation

% dual subgradient algorithm iterations
dualSubIters=2e3;

disp('generating problem ...')

% number of agents
AgN=3;  
% agent's dimension
Agni=4;

% problem data generation, see progen.m
probdata=progen(AgN,Agni);

% graph characteristics
stocasticity='doubly';
graphtype='binomial';
graphprob=1;% €(0,1], for binomial only
% graph generation, see myGraph.m
[ANW,A]=myGraph(AgN,stocasticity,graphtype,graphprob);
 
disp('generated!')

% problem data assignation
b=probdata.b;c=probdata.c;
d=probdata.d;D=probdata.D;
H=probdata.H;UB=probdata.UB;
LB=probdata.LB;

%% dual subgradient resolution
disp('solving with dual subgradient method...')

[primal_cost,dual_cost,primal_cost_RA,dual_cost_RA,consensus_error,...
    lambda,ZZ,ZRA] =dualsub(AgN,Agni,A,ANW,b,c,d,D,H,LB,UB,dualSubIters);
%d,D,LB unused

%% centralized solution
disp('centralized solution');

% linprog options
options = optimoptions('linprog','Display','none');

[~, centr_cost, exit_flag] = linprog(c,D,d,H,b,LB,UB,options);

% error catching
if exit_flag ~= 1
  fprintf(2,'A problem occurred in the centralized solution\n');
  return;
end

fprintf('Centralized optimal cost is %.4g\n',centr_cost);

%% results
disp('printing...')

% cost plot
figure
  % semilog plot for y axis for primal cost error
  semilogy(1:dualSubIters,abs(primal_cost(1:dualSubIters)-centr_cost), 'LineWidth',2);
  hold on, grid on, zoom on
  % for running average of primal cost  error
  semilogy(1:dualSubIters,abs(primal_cost_RA(1:dualSubIters)-centr_cost), 'LineWidth',2);
  % for dual cost  error
  semilogy(1:dualSubIters,abs(dual_cost(1:dualSubIters)-centr_cost), 'LineWidth',3.5);
  % for running average of dual cost  error
  semilogy(1:dualSubIters,abs(dual_cost_RA(1:dualSubIters)-centr_cost), 'LineWidth',2);
  xlabel('t')
  ylabel('cost error')
  legend('primal cost','primal cost with running avg',...
                'dual cost','dual cost with running avg')

% 
zum=reshape(sum(ZZ),Agni,dualSubIters);
zumra=reshape(sum(ZRA),Agni,dualSubIters);
for kk=1:Agni
    figure
      plot(1:dualSubIters,zum(kk,:)-b(kk), 'LineWidth',2);
      hold on, grid on, zoom on
      plot(1:dualSubIters,zumra(kk,:)-b(kk), 'LineWidth',2);
      xlabel('t')
      ylabel('x_1^t + ... + x_N^t - b')
      legend('z','z from running avg')
end  

% consensus error, multiple nested plots
consMatr4print=reshape(consensus_error,Agni,dualSubIters);
figure;grid on
% tiledlayout(Agni,1);%% method1
for kk=1:Agni
    consRow4print=consMatr4print(kk,:);
%     nexttile%% method1
    subplot(Agni,1,kk);%% method2
    plot(1:dualSubIters,consRow4print(1:dualSubIters), 'LineWidth',2);
    grid on, zoom on
end
sgtitle('consensus error')