%% monte carlo results
% 
% DISTRIBUTED CONTROL SYSTEMS course
% 
% GROUP 22
% CANELLO GIANMARCO
% CERRI FRANCESCO
% RONCATO MARCO
% 
% this code take as input some data from the monte carlo solver program and
% shows statistics on convergence and consensus errors

function [Lcentr_cost,Ldual_cost,Ldual_cost_RA,Lprimal_cost,...
            Lprimal_cost_RA,LconsError_avg,...
           LprimRA_relErr,LdualRA_relErr...
             ]=mcShow(dataFolder,pathMode)

this_scriptName='mcShow';%mfilename
%mfilename works only when running script from matlab command window

% add this script's location directory to path
scriptPath=fileparts(which(this_scriptName));
addpath(scriptPath);
        
% check dataFolder parameter
if isfolder(dataFolder)
    
    % dataFolder given wrt this script path
    if nargin<2 || strcmp(pathMode,'relative')
        dataPath=scriptPath+"/"+dataFolder;
    % dataFolder given as absolute path
    elseif strcmp(showMode,'absolute')
        dataPath=dataFolder;
    else
        warning("%s not understood; only absolute or relative allowed",pathMode);
    end
    
    % add dataPath to execution path
    addpath(dataPath);
    % get the list of .mat data element in dataFolder
    matFilesList=dir(dataFolder+"/*.mat");    
else
    error("%s is not a valid folder",dataFolder);
end


%% data aggregation

% mean values on repetition of same experiment

% retrieving number of repetitions
if contains(dataPath,"rep")
    repNumber=char(dataPath);
    repNumber=str2double(repNumber(strfind(dataPath,"rep")+3));
elseif contains(dataPath,"REP")
    repNumber=char(dataPath);
    repNumber=str2double(repNumber(strfind(dataPath,"REP")+3));
elseif contains(dataPath,"Rep")
    repNumber=char(dataPath);
    repNumber=str2double(repNumber(strfind(dataPath,"Rep")+3));  
else
    repNumber=1;
end

repNumber=10;

% retrieving number of dual subgradient iterations from a data field
sampleLen=length(matFilesList);
loadedMatData=load(matFilesList(end).name);
arrayLen=size(loadedMatData.consensus_error,3);
clear loadedMatData

% init of arrays based on dimension retrieved
Lcentr_cost=zeros(sampleLen/repNumber,1);% centralized cost
Ldual_cost=zeros(sampleLen/repNumber,arrayLen);% dual cost
Ldual_cost_RA=zeros(sampleLen/repNumber,arrayLen);% dual cost running average
LdualRA_relErr=zeros(sampleLen/repNumber,arrayLen);%d.c. r.a. relative error
Lprimal_cost=zeros(sampleLen/repNumber,arrayLen);% primal cost
Lprimal_cost_RA=zeros(sampleLen/repNumber,arrayLen);%primal c. running average
LprimRA_relErr=zeros(sampleLen/repNumber,arrayLen);% p.c. r.a. relative error
LconsError_avg=zeros(sampleLen/repNumber,arrayLen);% averaged consensus error
%NOTE reduced tensor dimension to a bidimensional array, averaging on the
%"lambdas'" dimension

%reduced name list parameters, excluding repetitions "REPn"
newNameList="";
mypat="";


% data averaging on multiple, similar, experiments
for loadCounter=1:length(matFilesList)
    % counter for averaging considering the number of repeated experiments
    reduCounter=fix((loadCounter+repNumber-1)/repNumber);
    % running average counter
    runCounter=1+rem((loadCounter+repNumber-1),repNumber);
    % load data fields from .mat files in the given folder
    loadedMatData=load(matFilesList(loadCounter).name);
        % centarized
        Lcentr_cost(reduCounter)=myRunAve(Lcentr_cost(reduCounter,:),loadedMatData.centr_cost,runCounter);
        
%         % dual
        Ldual_cost(reduCounter,:)=myRunAve(Ldual_cost(reduCounter,:),(loadedMatData.dual_cost)',runCounter);
        
        Ldual_cost_RA(reduCounter,:)=myRunAve(Ldual_cost_RA(reduCounter,:),(loadedMatData.dual_cost_RA)',runCounter);
        
        LdualRA_relErr(reduCounter,:)=myRunAve(LdualRA_relErr(reduCounter,:),...
            myRelErr((loadedMatData.dual_cost_RA)',Lcentr_cost(reduCounter)),runCounter);

%         % primal
        Lprimal_cost(reduCounter,:)=myRunAve(Lprimal_cost(reduCounter,:),(loadedMatData.primal_cost)',runCounter);
        
        Lprimal_cost_RA(reduCounter,:)=myRunAve(Lprimal_cost_RA(reduCounter,:),(loadedMatData.primal_cost_RA)',runCounter);
        
        LprimRA_relErr(reduCounter,:)=myRunAve(LprimRA_relErr(reduCounter,:),...
            myRelErr((loadedMatData.primal_cost_RA)',Lcentr_cost(reduCounter)),runCounter);

%         % consensus error
        LconsError_avg(reduCounter,:)=myRunAve(...
            LconsError_avg(reduCounter,:),reshape(mean(loadedMatData.consensus_error),[1,arrayLen]),runCounter);

        % reduced name list extraction, leaving out repeated experiments names
        newName=matFilesList(loadCounter).name(1:end-4);
        if strcmp(mypat,"")% init, first element to search
            newNameList(1,1)=newName;%+".mat";
            mypat=newName;% pattern to search
            fillCount=2;% filling counter
        elseif not(contains(newName,mypat))% new name found
            newNameList(fillCount,1)=newName;%+".mat";
            mypat=newName;
            fillCount=fillCount+1;
        end
end

clear fillCount mypat newName  runCounter% reduCounter

%% data aggregation 
% data of experiments with a common characteristic (number of agents N or
% agents' dimension M) represented together



% based on costant agents' dimension
dimQuery="9";
searchQuery="M"+dimQuery;%will search experiments with

% figures for plotting dual/primal, running averages costs and consensus
% error
mc_primDual_M_plot=figure(1);
% set(mc_primDual_M_plot,'Visible','off');
legarray1="";% array for constructing the legend
mc_cons_M_plot=figure(2);
% set(mc_cons_M_plot,'Visible','off');
legarray2="";
legCounter=1;%dimension of legend
grid on, zoom on

warning('off','MATLAB:handle_graphics:exceptions:SceneNode');

for queryCount=1:reduCounter
    % searching the requested characteristic
    if contains(newNameList(queryCount),searchQuery)
%         fprintf("%s\n",newNameList(queryCount));
        % legend construction
        legName=char(newNameList(queryCount));
        
        legarray1(legCounter:legCounter+1)=...
          ["dual "+legName(1:2)+","+legName(4:6) "primal "+legName(1:2)+","+legName(4:6)];
        % plot
        figure(1);hold on, zoom on, grid on
%         semilogy(1:arrayLen,abs(Ldual_cost_RA(queryCount,1:arrayLen)-Lcentr_cost(queryCount)), 'LineWidth',2);
%         semilogy(1:arrayLen,abs(Lprimal_cost_RA(queryCount,1:arrayLen)-Lcentr_cost(queryCount)), 'LineWidth',2);
        semilogy(1:arrayLen,LdualRA_relErr(queryCount,1:arrayLen), 'LineWidth',2);
        semilogy(1:arrayLen,LprimRA_relErr(queryCount,1:arrayLen), 'LineWidth',2);
        sgtitle('error of dual and primal costs by running average');xlabel('t')
        hold off
        
        % same as above, for consensus plot
        figure(2);hold on, zoom on, grid on        
        legarray2(legCounter:legCounter+1)=(legName(1:2)+","+legName(4:6));
        legCounter=legCounter+2;
        semilogy(1:arrayLen,LconsError_avg(queryCount,1:arrayLen) ,'LineWidth',2);
        sgtitle('consensus error');xlabel('t')
        hold off        
    end

end
% legend assignation
warning('off','MATLAB:legend:IgnoringExtraEntries');
figure(1);legend(legarray1);
figure(2);legend(legarray2);
warning('on','MATLAB:legend:IgnoringExtraEntries');        
        


% based on costant number of agents
% same as above
agenQuery="5";
searchQuery="N"+agenQuery;% will search for experiments with

mc_primDual_N_plot=figure(3);
% set(mc_primDual_N_plot,'Visible','off');
legarray3="";
mc_cons_N_plot=figure(4);
% set(mc_cons_N_plot,'Visible','off');
legarray4="";
legCounter=1;
grid on, zoom on

for queryCount=1:reduCounter
    if contains(newNameList(queryCount),searchQuery)
%         fprintf("%s\n",newNameList(queryCount));
        legName=char(newNameList(queryCount));
        
        legarray3(legCounter:legCounter+1)=...
          ["dual "+legName(1:2)+","+legName(4:6) "primal "+legName(1:2)+","+legName(4:6)];
        figure(3);hold on, zoom on, grid on
%         semilogy(1:arrayLen,abs(Ldual_cost_RA(queryCount,1:arrayLen)-Lcentr_cost(queryCount)), 'LineWidth',2);
%         semilogy(1:arrayLen,abs(Lprimal_cost_RA(queryCount,1:arrayLen)-Lcentr_cost(queryCount)), 'LineWidth',2);
        semilogy(1:arrayLen,LdualRA_relErr(queryCount,1:arrayLen), 'LineWidth',2);
        semilogy(1:arrayLen,LprimRA_relErr(queryCount,1:arrayLen), 'LineWidth',2);
        sgtitle('error of dual and primal costs by running average');xlabel('t')
        hold off
        
        figure(4);hold on, zoom on, grid on        
        legarray4(legCounter:legCounter+1)=(legName(1:2)+","+legName(4:6));
        legCounter=legCounter+2;
        semilogy(1:arrayLen,LconsError_avg(queryCount,1:arrayLen) ,'LineWidth',2);
        sgtitle('consensus error');xlabel('t')
        hold off        
    end

end
warning('off','MATLAB:legend:IgnoringExtraEntries');
figure(3);legend(legarray3);    
figure(4);legend(legarray4);
warning('on','MATLAB:legend:IgnoringExtraEntries');
pause(1);% helps slower cpus avoid incurring in queued warnings
warning('on','MATLAB:handle_graphics:exceptions:SceneNode');


%% exporting data

% saving folder
saveFolder="MC_aggr";
warning('off','MATLAB:MKDIR:DirectoryExists');
mkdir(saveFolder);
warning('on','MATLAB:MKDIR:DirectoryExists');
% save filename and filepath
solution_filename="aggr_REP"+repNumber;
sol_filepath=saveFolder+"/"+solution_filename;

% exporting .mat file
save(sol_filepath+'_[M'+dimQuery+'_N'+agenQuery+']',"-regexp","Err");
save(sol_filepath+'_[M'+dimQuery+'_N'+agenQuery+']',"-append","-regexp","cost");
save(sol_filepath+'_[M'+dimQuery+'_N'+agenQuery+']',"-append","reduCounter");
save(sol_filepath+'_[M'+dimQuery+'_N'+agenQuery+']',"-append","arrayLen");
save(sol_filepath+'_[M'+dimQuery+'_N'+agenQuery+']',"-append","-regexp","new");
save(sol_filepath+'_[M'+dimQuery+'_N'+agenQuery+']',"-append","-regexp","repN");

% exporting plots as .jpg and .fig
%NOTE [term] not costants, but usefull fol file duplicates
saveas(mc_cons_M_plot,sol_filepath+'_ConsPlot_M'+dimQuery+'_[N'+agenQuery+']_.jpg');
saveas(mc_cons_M_plot,sol_filepath+'_ConsPlot_M'+dimQuery+'_[N'+agenQuery+']_.fig');
% close(mc_cons_M_plot);
% clear mc_cons_M_plot  

saveas(mc_primDual_M_plot,sol_filepath+'_CostPlot_M'+dimQuery+'_[N'+agenQuery+']_.jpg');
saveas(mc_primDual_M_plot,sol_filepath+'_CostPlot_M'+dimQuery+'_[N'+agenQuery+']_.fig');
% close(mc_primDual_M_plot);
% clear mc_primDual_M_plot  

saveas(mc_primDual_N_plot,sol_filepath+'_CostPlot_N'+agenQuery+'_[M'+dimQuery+']_.jpg');
saveas(mc_primDual_N_plot,sol_filepath+'_CostPlot_N'+agenQuery+'_[M'+dimQuery+']_.fig');
% close(mc_primDual_N_plot);
% clear mc_primDual_N_plot    

saveas(mc_cons_N_plot,sol_filepath+'_ConsPlot_N'+agenQuery+'_[M'+dimQuery+']_.jpg');
saveas(mc_cons_N_plot,sol_filepath+'_ConsPlot_N'+agenQuery+'_[M'+dimQuery+']_.fig');
% close(mc_cons_N_plot);
% clear mc_cons_N_plot 

end



% relative error 
function relErr=myRelErr(vect_rel,vect_base)
    relErr=abs(vect_rel-vect_base)/vect_base;
end

% running average
function runAve=myRunAve(ave2run,newInput,runCounter)
    runAve=ave2run*(runCounter-1)/runCounter+newInput/runCounter;
end
   