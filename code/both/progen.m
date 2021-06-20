function progenres=progen(N,M,toggleTask2,areaWidth,...
                    areaHeigth,noiseOn,noRestrictions,spawnSeed)

%% problem generator for constrained coupled linear problems
% 
%  TASK1:
%   input: number of agents N, agent's dimension M (M>N)
%NOTE: generator code is imported from the one by A.Testa

%  TASK2:
%   input: number of agents/tasks N, M (M=N), Task2 mode activation toggleTask2,
%  spawn area width and height areaWidth areaHeight, toggle noise on cost noiseOn,
%  no restrictions on which agent can perform which task, rng seed spaqnSeed
%  
%   output: see each task
%% generation
% task 2
%    at least 2 arguments passed and toggled 
    if nargin>2 && toggleTask2
        
% 	 same number of tasks and agents
        M=N;

%        check if spawnSeed vector has enough elements, otherwise fills it
        if length(spawnSeed)<2,spawnSeed=randi(10,[1,2]); end
        
% 	     randomly(seeded) generation of agents and tasks on a given area
        agents=spawnEntities(false,N,areaWidth,areaHeigth,spawnSeed(1));
        tasks=spawnEntities(false,N,areaWidth,areaHeigth,spawnSeed(2));

        %     ci= ci (tasks', robots' positions)+noise
%         cost noise
        noizCoeff=0.3;
%         noise exclusion
        if not(noiseOn), noizCoeff=0; end
        
% 	 cost for an agents to perform a task, based on the distance in between
        c=zeros(N*N,1);
        for ii=1:N
            for kk=1:N
                rnoiz=noizCoeff*rand();

                c((ii-1)*N+kk)=sqrt((agents(ii,1)-tasks(ii,1))^2+...
                                 (agents(ii,2)-tasks(ii,2))^2)+rnoiz;
            end
        end
        

%   task assignment matrix, N repeated (for each agent) M*M (for each task)
%     Hi subass(eye(n)) | missing column k if agent i cant perform task k
%           Hi€R^NxKi
%           Hi=Hi(tasks, robots):1) evaluation of distances (>max)
%                                2) robot type that can perform a task
% 	 if true, toggles off limitations on which agent may do each task
        if noRestrictions
            H=repmat(eye(M), [1,N]); 
        else
            H=zeros(N,N*M);
%            tries to randomly create a matrix that has each task servable
%             by at least an agent
            while 1
                for ii=1:N
                    Hi=diag(round(1.49*rand(1,M)));
                    H(:,(ii-1)*N+1:(ii-1)*N+M)=Hi;
                end
                if sum(H,2)>=1, break,end%servability test
            end
        end

% constraint on assignation: each task should be served by a single agent
% and viceversa
        b = ones(M,1);

% constraints
        G=zeros(N,N^2);
        for ii=1:N
            G(ii,(ii-1)*N+1:ii*N)=ones(1,N);
        end
        g=ones(N,1);

        
%     Di=0
%     di=0

% 	 bounds for agents' states
        LB=zeros(M*N,1);
        UB=ones(M*N,1);
%  structured data
        progenres.c=c; progenres.b=b;
        progenres.G=G;progenres.g=g;
        progenres.H=H;
        progenres.LB=LB;progenres.UB=UB;
        progenres.agents=agents; progenres.tasks=tasks;
        progenres.areaWidth=areaWidth;
        progenres.areaHeigth=areaHeigth;
        
%#####################################################################
% selects task 1

   else
    % N agents of size M, M>N
    % ci € R^ni
    % Hi € R^S x ni
    % b € R^S
    % Pi = {zi € R^ni | Di*zi=<di, Gi*zi=gi}
    % Di € R^ ni x ni, di € R^ni
    % Gi € R^ ni x ni, gi € R^ni
    % qui:
    % ci € R^M, c € R^M*N
    % Hi € R^M x M , H € R^M x M*N
    % b € R^M
    % Pi | Gi=0, Di € R^N x N*M
    % -->
    % ni=M
    % S=M
   
        if not(M>N),error('mandatory: M>N, N first value'),end

% 	 profit associated with assigning job j to machine i
        P_ = randi([10 25],N,M);
        P  = max(P_,[],'all') +1 - P_ ;

% 	 claim on the capacity of machine i by job j
        W = randi([5 25],N,M);
        d = zeros(N,1); % d=[...dj...]

        for ijk = 1:N 
            J_i = [];
            for jjj = 1:M
              p_j = P(:,jjj);
              [~, ind_m] = min(p_j);
              ind_ = find(ind_m == ijk);
              if ~isempty(ind_)
                J_i = [J_i jjj];
              end
            end
            c_v = zeros(N,1);
            if ~isempty(J_i)
                for lll = 1:N
                    w_l = W(lll,:);
                    for jkl=1:numel(J_i)
                        c_v(lll)=c_v(lll)+w_l(J_i(jkl));
                    end
                end
            end
            d(ijk) = 9*(M/N) + 0.4*max(c_v); 
        end

        b = ones(M,1); %M=SS
        H = repmat(eye(M,M), [1,N]);
       
%	 machine capacity, local constraints
        D = zeros(N,M*N);% D=[...Dj...]
        for i = 1:N
            D(i,M*(i-1)+1:M*(i-1)+M) = W(i,:);
        end
         
%	 maximize profit
        c = reshape(P',[1,M*N]);% c=[...cj...]
        LB = zeros((N*M),1);
        UB = ones((N*M),1);


%       structured data
        progenres.b=b;progenres.c=c;progenres.d=d;progenres.D=D;
        progenres.H=H;progenres.LB=LB;progenres.UB=UB;

    end
end
