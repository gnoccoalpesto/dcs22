%% monte carlo experiments for assignatino task
% 
% DISTRIBUTED CONTROL SYSTEMS course
% 
% GROUP 22
% CANELLO GIANMARCO
% CERRI FRANCESCO
% RONCATO MARCO
% 
% this code take a range of conditions for [...]
% generates the data and solves the problem using dual subgradient distributed
% algorithm. 
% results are outputted in form of .mat, .txt
% plots in form of .fig, .jpg
close all
clear
clc
rng(69);

this_scriptName='task2MC';%mfilename
%mfilename works only when running script in matlab command window
    
addpath(fileparts(which(this_scriptName)));
addpath(fullfile(fileparts(which(this_scriptName)),'/both'));
addpath(fullfile(fileparts(which(this_scriptName)),'/task2'));
    

%% data range

% dual subgradient algorithm iterations
dualsubIter=8e3;
% number of agents
AgN_lib=3:14;
% agents dimension: same as number of agents
% stocasticity
stocasticity_lib="doubly";%["row","doubly"];
%graph topology
graphTopology_lib="binomial";%["binomial","path","cyclic"];
% graph probability (binomial only)
graphProbability_lib=1;%0.1*(5:10);
% spawn area characteristichs
areaWidth=1;
areaHeight=1;
% repetitions of all monte carlo experiment
numRedo=8;
% rngseeds
agentSeed=0;
taskSeed=0;
%NOTE 0 == use time

%total number of experiments
experTot=length(AgN_lib)*length(stocasticity_lib)*numRedo...
    (length(graphTopology_lib)-1+length(graphProbability_lib));
    
% current experiment counter
experCount=0;
% failed experiments counter
warningCount=0;

% linprog options (no message displayed)
linprog_options = optimoptions('linprog','Display','none');

% results' data folder
data_mainFolder="MCassignment";
mkdir(data_mainFolder);
% %subfolders
% figFolder
% txtFolder
% matFolder
% jpgFolder

% path to log
logFile=data_mainFolder+"/log_monteCarlo.txt";

% iteration times
% iterTime=zeros(experTot);

%% monte carlo experiments

% experiments iteration
for AgN=AgN_lib
for stocasticity=stocasticity_lib
for graphTopology=graphTopology_lib
    nobin_iter_count=0;%execution counter (see below)
for graphProbability=graphProbability_lib
% executes once or if binomial graph
if strcmp(graphTopology,"binomial") || nobin_iter_count<1
for iterRedo=1:numRedo   
    % iteration progression (*/tot numb. iterations)
    experCount=experCount+1;
    % execution counter+=1
    nobin_iter_count=nobin_iter_count+1;
    
    % experiment characteristics
    iterSet_r1="number of agents: "+AgN;
    iterSet_r2="stocasticity: "+stocasticity+", graph topology: "+graphTopology;
    if strcmp("binomial",graphTopology)
        iterSet_r3="binomial graph probability: "+graphProbability+"\n";
    else
        iterSet_r3="";
    end
    if iterRedo>1
        iterSet_r4="repeted "+(iterRedo-1)+" times after first iteration      \n";
    else
        iterSet_r4="";
    end
    
    fprintf("iteration:  %d/%d\n"+iterSet_r1+"\n"+iterSet_r2+"\n"+iterSet_r3+iterSet_r4+"\n",experCount,experTot);

    % nominal execution
    try
        % data generation
        Agni=AgN;
        
        msglen1=fprintf("generating problem ...");
        probdata=...
            progen(AgN,Agni,true,areaWidth,areaHeight,true,true,[agentSeed,taskSeed]);
        % graph
        [ANW,A]=myGraph(AgN,stocasticity,graphTopology,graphProbability);
        % problem
        G=probdata.G;H=probdata.H;b=probdata.b;
        c=probdata.c;g=probdata.g;
        UB=probdata.UB;LB=probdata.LB;
        areaHeigth=probdata.areaHeigth;
        areaWidth=probdata.areaWidth;
        Ag=probdata.agents;
        Ts=probdata.tasks;
        
        msglen1=msglen1+fprintf(" done!\n");
        
        % distributed solution
        msglen1=msglen1+fprintf("solving with dual subgradient method...\n");
        [primal_cost,dual_cost,primal_cost_RA,dual_cost_RA, consensus_error,...
           ~,ZZ,ZRA]= twodualsub(dualsubIter,AgN,A,ANW,b,c,g,G,H,LB,UB);
        
        % centralized solution
        [~, centr_cost, exit_flag] = linprog(c,[],[],[H;G],[b;g],LB,UB,linprog_options);
        
        % assignment        
        Ass_mat = round(ZRA(:,:,dualsubIter))>0.85;
        [ag2assign,ts4assign]=find(Ass_mat==1);
        
        msglen1=msglen1+fprintf("saving data... ");
        
        % savefile path creation using experiment's characteristichs
        solution_filename="assign_N"+AgN+"seed_["+agentSeed+"_"+taskSeed+"]";
        if iterRedo>1
            solution_filename=solution_filename+"_rep"+(iterRedo-1);
        end
        sol_filepath=data_mainFolder+"/"+solution_filename;
        
        % exporting data to mat file
        save(sol_filepath,"-regexp","cost");
        save(sol_filepath,"-append","consensus_error");
        save(sol_filepath,"-regexp","ZZ","-append");
        save(sol_filepath,"-struct","probdata","-append");
        save(sol_filepath,"-regexp","ass","-append");
        save(sol_filepath,"Ag","-append");
        save(sol_filepath,"Ts","-append");
                
        % exporting cost difference plot as jpg and fig
        %NOTE when imported .fig data has still visibility off; must use
        %command: set(FIGNAME,'Visible','on') and will show up on the spot
        
        % assignation plot
        assiPlot=figure();set(assiPlot,'Visible','off');
        plot(Ag(:,1),Ag(:,2),'go');
        hold on
        pause(1);
        plot(Ts(:,1),Ts(:,2),'rx')
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
        sgtitle('assignation plot')
        hold off
        saveas(assiPlot,sol_filepath+'_assiPlot.jpg');
        saveas(assiPlot,sol_filepath+'_assiPlot.fig');
        close(assiPlot);
        clear assiPlot 
        % assignment cost error plot
        
        % cost plot
        primal_dual_plot=figure();
        set(primal_dual_plot,'Visible','off');
        semilogy(1:dualsubIter,abs(primal_cost(1:dualsubIter)-centr_cost), 'LineWidth',2);
        hold on, grid on, zoom on
        semilogy(1:dualsubIter,abs(primal_cost_RA(1:dualsubIter)-centr_cost), 'LineWidth',2);
        semilogy(1:dualsubIter,abs(dual_cost(1:dualsubIter)-centr_cost), 'LineWidth',4);
        semilogy(1:dualsubIter,abs(dual_cost_RA(1:dualsubIter)-centr_cost), 'LineWidth',2);
        xlabel('t')
        ylabel('cost error')
        legend('primal cost','primal cost with running avg',...
                    'dual cost','dual cost with running avg','Location','northoutside')
        hold off
        saveas(primal_dual_plot,sol_filepath+'_CostPlot.jpg');
        saveas(primal_dual_plot,sol_filepath+'_CostPlot.fig');
        close(primal_dual_plot);
        clear primal_dual_plot   
        
        
%         % consensus error plot
%         consMatr4print=reshape(consensus_error,Agni,dualsubIter);
%         consErr_plot=figure();
%         set(consErr_plot,'Visible','off');
%         for kk=1:Agni
%             consRow4print=consMatr4print(kk,:);
%             subplot(Agni,1,kk);
%             plot(1:dualsubIter,consRow4print(1:dualsubIter), 'LineWidth',2);
%             grid on, zoom on
%         end
%         sgtitle('consensus error')
%         saveas(consErr_plot,sol_filepath+'_consensusPlot.jpg');
%         saveas(consErr_plot,sol_filepath+'_consensusPlot.fig');
%         close(consErr_plot);
%         clear consErr_plot 
        
        % clean up for next iteration
        msglen1=msglen1+fprintf("done!");
        
        
%             clear -regexp cost probdata
%             clear -regexp 
        % this will delete some command window lines
        fprintf((repmat('\b', 1,msglen1+1)));
            
    % error catching in case of failed experiment
    catch
        warning("solution failed. Exception written in log file: %s\n",logFile);
        % log message
        logMsg=iterSet_r1+iterSet_r2;
        if strcmp("binomial",graphTopology)
            logMsg=logMsg+iterSet_r3;
        end
        % logger
        writematrix('experiment',logFile,'WriteMode','append');
        writematrix(logMsg,logFile,'WriteMode','append');
        writematrix('at time',logFile,'WriteMode','append');
        writematrix(clock,logFile,'WriteMode','append');
        writematrix('################',logFile,'WriteMode','append');
        warningCount=warningCount+1;
    end
%     iterTime(experCount)=toc;
%     fprintf("elapsed time:\n-this iteration: %f\n--since start: %f\n",...
%     iterTime(experCount),sum(iterTime));
    disp('###########');
end
end
end
end
end
end
close all

% results show off
% fprintf("total elapsed time: %f\naverage time for experiment: %f\n"...
%     +"min: %f\nmax: %f\n",sum(iterTime),mean(iterTime),min(iterTime),max(iterTime));
if warningCount==0
    fprintf("_-^-_ all the experiments ended successfully! _-^-_\n");
else
    warning("a total of %d experiments failed. Check %s\n",warningCount,logFile);
end