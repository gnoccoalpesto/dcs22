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
AgN=3;
%dimension of agents
Agni=4;

%spawn agenti

% spawn task

% generazione problema 2
% c=c(distanza agente i task k)=cik  C€R^N*Ki:sottoassieme(R^NxN)
for ii=1:N %agents
    %Ki=N_In(k);
    for kk=1:N%==Ki
%         cik=distfunc(sqrt((xagent-xtask)^2+(yagent-ytask)^2));
        cik=(sqrt((xagent-xtask)^2+(yagent-ytask)^2));
%         cik=cik+random_noise;
    end
end
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