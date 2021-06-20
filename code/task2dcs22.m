%% DCS 2020
% PROJECT 3
% GROUP 22: CANELLO
%           CERRI
%           RONCATO

%TASK 2
% generation of an assignment problem of N tasks to N agents and resolution
% by the means of a distributed dual subgradient algorithm

clear,clc,close all
rng(69);

% adds to path the folder containing this file and library subfolders
addpath(fileparts(which('task2dcs22')));
addpath(fullfile(fileparts(which('task2dcs22')),'/both'));
addpath(fullfile(fileparts(which('task2dcs22')),'/task2'));


%% problem generation

disp('generating problem ...')

% number of agents
AgN=8;
% dimension of agents
Agni=AgN; % == number of tasks == agent number
% number of iterations
maxIters=5e1;

% graph characteristics
stocasticity='doubly';
% graph topology
graphtype='binomial';
graphprob=1;% €(0,1]
% graph generation
[ANW,A]=myGraph(AgN,stocasticity,graphtype,graphprob);

% problem data generation
probdata=progen(AgN,Agni,true,1,1,true,true,[5,3]);
% data destructuration
G=probdata.G;H=probdata.H;b=probdata.b;
c=probdata.c;g=probdata.g;
UB=probdata.UB;LB=probdata.LB;
areaHeigth=probdata.areaHeigth;
areaWidth=probdata.areaWidth;
Ag=probdata.agents;
Ts=probdata.tasks;

disp('generated!')

%% distributed dual subgradient method solution
disp('solving with dual subgradient method...    ')

[primal_cost,dual_cost,...
 primal_cost_RA,dual_cost_RA,...
 consensus_error,lambda,ZZ,ZRA   ]...
                        = twodualsub(maxIters,AgN,A,ANW,b,c,g,G,H,LB,UB);

%% centralized solution

options = optimoptions('linprog','Display','none');
[~, centr_cost, exit_flag] = linprog(c,[],[],[H;G],[b;g],LB,UB,options);

%% results

disp('printing...')

% plot of agents and tasks: tasks are red X marks, agents green circles
figure();
plot(Ag(:,1),Ag(:,2),'go');
hold on
pause(1);
plot(Ts(:,1),Ts(:,2),'rx');

% assignation thresholding
Ass_mat = round(ZRA(:,:,maxIters))>0.85;
%NOTE could have been done rounding to 1 too

% assignation extraction
[ag2assign,ts4assign]=find(Ass_mat==1);
% assignation display
  for ii=1:length(ag2assign)
      
      agent_x=Ag(ag2assign(ii),1);
      agent_y=Ag(ag2assign(ii),2);
      
      task_x=Ts(ts4assign(ii),1);
      task_y=Ts(ts4assign(ii),2);
      line([agent_x,task_x],[agent_y,task_y],'color',[rand rand rand]);
      
      
  end
% area boundaris reshaping (graphical purpuses only)
line([0,0],[0,areaHeigth],'color',[0 0 0]);
line([0,areaWidth],[0,0],'color',[0 0 0]);
line([areaWidth,areaWidth],[0,areaHeigth],'color',[0 0 0]);
line([0,areaWidth],[areaWidth,areaHeigth,],'color',[0 0 0]);
hold off

% assignment cost error plot
cost4assign=0;
for ii=1:AgN
    cost4assign=cost4assign+c((ag2assign(ii)-1)*AgN+ts4assign(ii));
end

figure()
semilogy(1:maxIters,abs(primal_cost(1:maxIters)-centr_cost), 'LineWidth',2);
hold on, grid on, zoom on
semilogy(1:maxIters,abs(primal_cost_RA(1:maxIters)-centr_cost), 'LineWidth',2);
semilogy(1:maxIters,abs(dual_cost(1:maxIters)-centr_cost), 'LineWidth',5);
semilogy(1:maxIters,abs(dual_cost_RA(1:maxIters)-centr_cost), 'LineWidth',2);
xlabel('t')
ylabel('cost error')
legend('primal cost','primal cost with running avg',...
                'dual cost','dual cost with running avg')


% consMatr4print=reshape(consensus_error,Agni,maxIters);
% for kk=11:Agni%works on single iteraction
%     figure()
%     consRow4print=consMatr4print(kk,:);
%       semilogy(1:maxIters,consRow4print(1:maxIters), 'LineWidth',2);
% %       hold on, grid on, zoom on
%       xlabel('t')
%       ylabel('consensus error')
% end