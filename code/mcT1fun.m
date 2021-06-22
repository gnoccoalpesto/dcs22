%% monte carlo solver
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

function mcT1fun(AgN_lib,Agni_lib,numRedo,dualsubIter,data_mainFolder)
this_scriptName='mcT1fun';%mfilename
%mfilename works only when running script in matlab command window
    
addpath(fileparts(which(this_scriptName)));
addpath(fullfile(fileparts(which(this_scriptName)),'/both'));
addpath(fullfile(fileparts(which(this_scriptName)),'/task1'));
    
% close all
% clear
% clc
%% data range

% dual subgradient algorithm iterations
% dualsubIter=5e3;
% number of agents
% AgN_lib=3:6;
% agents dimension: relative to AgN (constraint by problem formulation)
% Agnimin_diff=1;
% Agnimin_diff=Agni_diff(1);
% AgniMax_diff=Agni_diff(2);
% AgniMax_diff=3;
% stocasticity
stocasticity_lib="doubly";%["row","doubly"];
%graph topology
graphTopology_lib="binomial";%["binomial","path","cyclic"];
% graph probability (binomial only)
graphProbability_lib=1;%0.1*(5:10);


% repetitions of all monte carlo experiment
% numRedo=7;

%total number of experiments
experTot=length(AgN_lib)*length(stocasticity_lib)*numRedo*...
    (length(graphTopology_lib)-1+length(graphProbability_lib))...
       *length(Agni_lib);
    
% current experiment counter
experCount=0;
% failed experiments counter
warningCount=0;

% linprog options (no message displayed)
linprog_options = optimoptions('linprog','Display','none');

% results' data folder
if nargin<5
    data_mainFolder="MC_data";
end
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
for Agni=Agni_lib%AgN+(Agnimin_diff:AgniMax_diff)
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
    iterSet_r1="number of agents: "+AgN+", agents' dimension: "+Agni;
    iterSet_r2="stocasticity: "+stocasticity+", graph topology: "+graphTopology;
    if strcmp("binomial",graphTopology)
        iterSet_r3="binomial graph probability: "+graphProbability+"\n";
    else
        iterSet_r3="";
    end
    if iterRedo>1
        iterSet_r4="repeted "+(iterRedo-1)+" times after first iteration              \n";
    else
        iterSet_r4="";
    end
    
    fprintf("iteration:  %d/%d\n"+iterSet_r1+"\n"+iterSet_r2+"\n"+iterSet_r3+iterSet_r4+"\n",experCount,experTot);
    
%     tic;
    % nominal execution
    try
        % data generation
        msglen1=fprintf("generating problem ...");
        probdata=progen(AgN,Agni);
        % graph
        [ANW,A]=myGraph(AgN,stocasticity,graphTopology,graphProbability);
        % problem
        b=probdata.b;c=probdata.c;
        d=probdata.d;D=probdata.D;
        H=probdata.H;UB=probdata.UB;
        LB=probdata.LB;
        msglen1=msglen1+fprintf(" done!\n");
        
        % distributed solution
        msglen1=msglen1+fprintf("solving with dual subgradient method...\n");
        [primal_cost,dual_cost,primal_cost_RA,dual_cost_RA,consensus_error,...
            ~,ZZ,ZZRA] =dualsub(AgN,Agni,A,ANW,b,c,d,D,H,LB,UB,dualsubIter);
        % centralized solution
        [~, centr_cost, ~] = linprog(c,D,d,H,b,LB,UB,linprog_options);
        
        msglen1=msglen1+fprintf("saving data... ");
        
        % savefile path creation using experiment's characteristichs
        if strcmp(graphTopology,"binomial")
            solution_filename="N"+AgN+"_M"+Agni+"_STC"+stocasticity+...
             "_TOP"+graphTopology+"_PRB";
            if graphProbability==1                    
                solution_filename=solution_filename+graphProbability;
            else
                solution_filename=solution_filename+"0_"+10*graphProbability;
            end
        else
            solution_filename="N"+AgN+"_M"+Agni+"_STC"+stocasticity+...
             "_TOP"+graphTopology;
        end
        if iterRedo>1
            solution_filename=solution_filename+"_rep"+(iterRedo-1);
        end
        sol_filepath=data_mainFolder+"/"+solution_filename;
        
        % exporting data to mat file
        save(sol_filepath,"-regexp","cost");
        save(sol_filepath,"-append","consensus_error");
        save(sol_filepath,"-regexp","ZZ","-append");
        save(sol_filepath,"-struct","probdata","-append");
        
        % exporting data to txt file
        %NOTE .txt dimension may be to big to some readers, text rows
        %should be loaded progressively, not at once

%         writematrix('centralized cost',sol_filepath+".txt",'WriteMode','append');
%         writematrix(centr_cost,sol_filepath+".txt",'WriteMode','append');
%         writematrix('################',sol_filepath+".txt",'WriteMode','append');
%         writematrix("primal cost'",sol_filepath+".txt",'WriteMode','append');
%         writematrix(primal_cost',sol_filepath+".txt",'WriteMode','append');
%         writematrix('################',sol_filepath+".txt",'WriteMode','append');
%         writematrix("primal RA'",sol_filepath+".txt",'WriteMode','append');
%         writematrix(primal_cost_RA',sol_filepath+".txt",'WriteMode','append');
%         writematrix('################',sol_filepath+".txt",'WriteMode','append');
%         writematrix("dual cost'",sol_filepath+".txt",'WriteMode','append');
%         writematrix(dual_cost',sol_filepath+".txt",'WriteMode','append');
%         writematrix('################',sol_filepath+".txt",'WriteMode','append');
%         writematrix("dual RA'",sol_filepath+".txt",'WriteMode','append');
%         writematrix(dual_cost_RA',sol_filepath+".txt",'WriteMode','append');
%         writematrix('################',sol_filepath+".txt",'WriteMode','append');
%         writematrix('consensus error',sol_filepath+".txt",'WriteMode','append');
%         writematrix(consensus_error,sol_filepath+".txt",'WriteMode','append');
%         writematrix('################',sol_filepath+".txt",'WriteMode','append');
%         writematrix('ZZ',sol_filepath+".txt",'WriteMode','append');
%         writematrix(ZZ,sol_filepath+".txt",'WriteMode','append');
%         writematrix('################',sol_filepath+".txt",'WriteMode','append');
%         writematrix('ZZra',sol_filepath+".txt",'WriteMode','append');
%         writematrix(ZZRA,sol_filepath+".txt",'WriteMode','append');
%         writematrix('################',sol_filepath+".txt",'WriteMode','append');
%         writematrix("probdata: b'",sol_filepath+".txt",'WriteMode','append');
%         writematrix(probdata.b' ,sol_filepath+".txt",'WriteMode','append');
%         writematrix('probdata:c',sol_filepath+".txt",'WriteMode','append');
%         writematrix(probdata.c,sol_filepath+".txt",'WriteMode','append');
%         writematrix('probdata:d',sol_filepath+".txt",'WriteMode','append');
%         writematrix(probdata.d,sol_filepath+".txt",'WriteMode','append');
%         writematrix('probdata:D',sol_filepath+".txt",'WriteMode','append');
%         writematrix(probdata.D,sol_filepath+".txt",'WriteMode','append');
%         writematrix('probdata:H',sol_filepath+".txt",'WriteMode','append');
%         writematrix(probdata.H,sol_filepath+".txt",'WriteMode','append');
%         writematrix("probdata: upper bound'",sol_filepath+".txt",'WriteMode','append');
%         writematrix(probdata.UB',sol_filepath+".txt",'WriteMode','append');
%         writematrix("probdata: lower bound'",sol_filepath+".txt",'WriteMode','append');
%         writematrix(probdata.LB',sol_filepath+".txt",'WriteMode','append');
        
        
        % exporting cost difference plot as jpg and fig
        %NOTE when imported .fig data has still visibility off; must use
        %command: set(FIGNAME,'Visible','on') and will show up on the spot
        
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
        
        
        % consensus error plot
        consMatr4print=reshape(consensus_error,Agni,dualsubIter);
        consErr_plot=figure();
        for kk=1
            consRow4print=consMatr4print(kk,:);
            subplot(Agni,1,kk);
            plot(1:dualsubIter,consRow4print(1:dualsubIter), 'LineWidth',2);
            grid on, zoom on
        end
        sgtitle('consensus error')
        saveas(consErr_plot,sol_filepath+'_consensusPlot.jpg');
        saveas(consErr_plot,sol_filepath+'_consensusPlot.fig');
        close(consErr_plot);
        clear consErr_plot 
        
        % clean up for next iteration
        msglen1=msglen1+fprintf("done!");
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
    disp('###########');
end
end
end
end
end
end
end

% results show off
if warningCount==0
    fprintf("_-^-_ all the experiments ended successfully! _-^-_\n");
else
    warning("a total of %d experiments failed. Check %s\n",warningCount,logFile);
end

end