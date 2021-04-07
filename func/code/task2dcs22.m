clear,clc,close all
rng(1);

% DCS 2020
% PROJECT 3
% GROUP 22: CANELLO
%           CERRI
%           RONCATO
%task2


% clear;close all;clc;

% adds to path the folder containing this file and /func subfolder
% change the following if the filename changes
% addpath(fileparts(which('dcs20_group22')));
% addpath(fullfile(fileparts(...
%     which('dcs20_group22')),'/func'));
% on .m file mfilename can be used to retrieve it


%% 

disp('generating problem ...')

%number of agents
AgN=6;
%dimension of agents
Agni=Agn;

% c=c(distanza agente i task k)=cik  C€R^N*Ki:sottoassieme(R^NxN)

stocasticity='doubly';
graphtype='binomial';
graphprob=1;% €(0,1]
[ANW,A]=myGraph(AgN,stocasticity,graphtype,graphprob);
%ANW unused (since it's easily identifiable from A)


b=probdata.b;c=probdata.c;
g=probdata.g;G=probdata.G;
H=probdata.H;UB=probdata.UB;
LB=probdata.LB;

disp('generated!')
%% spawn agents and tasks (temporary)
agents


%%
disp('solving with dual subgradient method...')

[primal_cost,dual_cost,primal_cost_RA,dual_cost_RA,...consensus_error,
    lambda,ZZ,ZRA,maxIters] =dualsub(AgN,Agni,A,ANW,b,c,d,D,H,LB,UB);
%d,D,LB unused

%% confronto centralizato distribuito
% disp('centralized solution');
% 
% options = optimoptions('linprog','Display','none');
% [~, fopt, exit_flag] = linprog(c,D,d,H,b,LB,UB,options)





%% improv
% barrier func per evitare che i robot si compenitrino (robot vel costante,
%       linea retta, pivoting)

%
%

Ag=probdata.agents;
Ts=probdata.tasks;

figure();
plot(Ag(:,1),Ag(:,2),'go');
    hold on
    pause(1);
    plot(Ts(:,1),Ts(:,2),'rx');

% Ass_mat = ZRA(:,:,tt) ;
% maybe recast assign mat in bool

[ag4assign,ts4assign]=find(Ass_mat==1);
pause(1);
  for ii=1:NN
      
      agent_x=Ag(ag4assign(ii),1); %Assign(ii,1) agent to be selected
      agent_y=Ag(ag4assign(ii),2);
      
      task_x=Ts(ts4assign(ii),1); %Assign(ii,1) task to be assinged
      task_y=Ts(ts4assign(ii),2);
      line([agent_x,task_x],[agent_y,task_y],'color',[rand rand rand]);
      
      
  end
hold off
















