% DCS 2020
% PROJECT 3
% GROUP 22: CANELLO
%           CERRI
%           RONCATO

%task2
clear,clc,close all
rng(69);

% adds to path the folder containing this file and library subfolders
addpath(fileparts(which('task2dcs22')));
addpath(fullfile(fileparts(which('task2dcs22')),'/both'));
addpath(fullfile(fileparts(which('task2dcs22')),'/task2'));


%% 

disp('generating problem ...')

% number of agents
AgN=12;
% dimension of agents
Agni=AgN;
% number of iterations
maxIters=1e4;
% 
stocasticity='doubly';
% graph topology
graphtype='binomial';
graphprob=1;% €(0,1]

[ANW,A]=myGraph(AgN,stocasticity,graphtype,graphprob);

% data generation
probdata=progen(AgN,Agni,true,1,1,true,true,[3,5]);

G=probdata.G;H=probdata.H;b=probdata.b;
c=probdata.c;g=probdata.g;
UB=probdata.UB;LB=probdata.LB;
areaHeigth=probdata.areaHeigth;
areaWidth=probdata.areaWidth;
Ag=probdata.agents;
Ts=probdata.tasks;

disp('generated!')

%% distributed solution
disp('solving with dual subgradient method...')

[primal_cost,dual_cost,primal_cost_RA,dual_cost_RA,consensus_error,...
    lambda,ZZ,ZRA] = twodualsub(maxIters,AgN,A,ANW,b,c,g,G,H,LB,UB);

%% centralized solution

options = optimoptions('linprog','Display','none');
[~, fopt, exit_flag] = linprog(c,[],[],[H;G],[b;g],LB,UB,options);

%% results

disp('printing...')

figure();
plot(Ag(:,1),Ag(:,2),'go');
    hold on
    pause(1);
    plot(Ts(:,1),Ts(:,2),'rx');

% assignation thresholding
% Ass_mat = ZRA(:,:,maxIters)>=.7;
Ass_mat = round(ZRA(:,:,maxIters))>0.99999;

[ag2assign,ts4assign]=find(Ass_mat==1);

  for ii=1:length(ag2assign)
      
      agent_x=Ag(ag2assign(ii),1);
      agent_y=Ag(ag2assign(ii),2);
      
      task_x=Ts(ts4assign(ii),1);
      task_y=Ts(ts4assign(ii),2);
      line([agent_x,task_x],[agent_y,task_y],'color',[rand rand rand]);
      
      
  end
% area boundaris
line([0,0],[agent_y,task_y],'color',[0 0 0]);
line([0,areaHeigth],[agent_y,task_y],'color',[0 0 0]);
line([areaWidth,0],[agent_y,task_y],'color',[0 0 0]);
line([areaWidth,areaHeigth],[agent_y,task_y],'color',[0 0 0]);
hold off

% assignment cost
cost4assign=0;
for ii=1:AgN
    cost4assign=cost4assign+c((ag2assign(ii)-1)*AgN+ts4assign(ii));
end

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