%% distributed dual subgradient solver for constr.coupl. lin.probl.
%
%  DCS GROUP 22
%       CANELLO
%       CERRI
%       RONCATO

% input: number of agents NN, agents' dimension nni, weighted/unwei. graph
% matrices AA AANW, equality constraint bb, agent's cost cc, inequality 
% constraints dd DD, assignment constraint matrix HH, lower/upper bound LBB
% UBB, algorithm iterations MAXITERS

% output: primal cost primCost, dual cost dualCost, running averages cost 
% primRA dualRA, consensus error consErr, lagrange multipliers lam, agents'
% state and running average ZZ ZRA

function...
    [primCost,dualCost,primRA,dualRA,consErr,...
                lam,ZZ,ZRA]=...
                           dualsub(NN,nni,AA,AANW,bb,cc,dd,...
                                    DD,HH,LBB,UBB,MAXITERS)

    %linprog options
    lpoptions = optimoptions('linprog','Display','none');
    
    msglnt2=fprintf('iterations: %d  \n',MAXITERS);

    msglnt2=fprintf('initializing...')+msglnt2;
    % cost init
    primCost=zeros(MAXITERS,1); % cost € R
    dualCost=zeros(MAXITERS,1);
    primRA=zeros(MAXITERS,1);
    dualRA=zeros(MAXITERS,1);
    consErr=zeros(1,nni,MAXITERS);
    
    % agents init
    % agent zi € R^ni, zz=[...;zi;...] agent matrix € R^N x ni
    % ZZ=[zz_1st_iteraction, ... , zz_MAXITERS]
    ZZ=zeros(NN,nni,MAXITERS); 
    ZRA=zeros(NN,nni,MAXITERS);
    
    SS=nni;
    % when %==SS is displayed, means that it's our particular case nni==SS
    %NOTE that this is a particular, easier case
    
    %multiplier init
    lam=zeros(NN,SS,MAXITERS);% 
    vv=zeros(NN,SS); % vv=[...;vi;...], vi€ R ^ SS
    
    % b splitting
    bi=bb/NN;
    %NOTE that this is a particular case; other eg: bi=1 if i=1, 0 if i!=1
    
    msglnt2=fprintf(" done!\n")+msglnt2;
    msglnt1=fprintf('progress: ');
    msglnt3=fprintf("0%%");
    
    % algorithm loop
    for tt = 1:MAXITERS-1
        % algorithm progression printing
        if mod(tt,MAXITERS/100)==0
            fprintf((repmat('\b', 1,msglnt3)));
            msglnt3=fprintf('%d%%',tt*100/MAXITERS);
        end
        
        %stepsize
        alpha_t = 1.6*(1/tt)^0.69; % 1/tt^alpha with alpha in (0.5, 1]
        % other examples of possible stepsizes
%         aa=0.1;%bb=0.1;% a,b>0
%         alpha_t=aa/(bb+tt);%square summable, not summable s.size
%         alpha_t=aa/sqrt(tt);%nonsummable diminishing s.size
%         alpha_t=aa/modulus(g(k),2);%nons. dim. const. s.length
                   

        for ii=1:NN
            % finding neighboor
            N_ii = find(AANW(:,ii) == 1)';            
            % update: averaging
            for kk=1:SS 
                vv(ii,kk) = AA(ii,ii)*lam(ii,kk,tt);
                % neighboor communication
                for jj = N_ii
                    vv(ii,kk) = vv(ii,kk) + AA(ii,jj)*lam(jj,kk,tt);
                end
            end
        end

        for ii=1:NN
            % data extraction for the current agent, based on dimension
            ci=cc((ii-1)*nni+1 : ii*nni);
            Hi=HH(:,nni*(ii-1)+1:nni*ii);
            UBi=UBB((ii-1)*nni+1:ii*nni);
            LBi=LBB((ii-1)*nni+1:ii*nni);
            Di=DD(ii,SS*(ii-1)+1:SS*ii);
            

            % minimization of agent quantity (solution of local problem)
            [ZZ(ii,:,tt),~,~]=linprog(ci+(vv(ii,:))*Hi', Di,dd(ii),...
            [],[],LBi,UBi,lpoptions);

            % Running average
            if tt==1
                ZRA(ii,:,tt) = ZZ(ii,:,tt);
            else
                ZRA(ii,:,tt) = (1/tt)*((tt-1)*ZRA(ii,:,tt-1)...
                                                  +ZZ(ii,:,tt));
            end
        end
    
    % multiplier update
        for ii=1:NN
            
            Hi=HH(:,nni*(ii-1)+1:nni*ii);
            % gradient computation
            grad_ii=Hi*ZZ(ii,:,tt)'-bi;
            %gather la_j^t from neigh | j€N_i
            lam(ii,:,tt+1)=vv(ii,:)'+alpha_t*grad_ii;
        end
        
        % multiplier average
        lamAvg = mean(lam(:,:,tt));

        % cost and consensus error computation
        for ii=1:NN
            
            ci=cc((ii-1)*nni+1 : ii*nni);
            
            % agents' primal cost contribute and sum
            ffii = ci * ZZ(ii,:,tt)';
            primCost(tt) = primCost(tt) + ffii;
            %
            ffii = ci * ZRA(ii,:,tt)';
            primRA(tt) = primRA(tt) + ffii;
            
            % agents' dual cost contribute and sum
            qqi = ci * ZZ(ii,:,tt)'...
                +lam(ii,:,tt) * ( ZZ(ii,:,tt) - bi' )';
            dualCost(tt) = dualCost(tt) + qqi;
            %
            qqi = ci *ZRA(ii,:,tt)'...
                + lam(ii,:,tt) * ( ZRA(ii,:,tt) - bi' )';
            dualRA(tt) = dualRA(tt) + qqi;

            % consensus error
            consErr(:,:,tt) = consErr(:,:,tt) + abs(lam(ii,:,tt) - lamAvg);
        end   
    end
    
    % last iteration
    tt=MAXITERS;
    
    % algorithm progress printing
    fprintf((repmat('\b', 1,msglnt1+msglnt3)));%display cleaning
    msglnt1=fprintf('progress: 100%%\n');
    
    % update averaging
    for ii=1:NN
        for kk=1:SS 
            
            N_ii = find(AANW(:,ii) == 1)';
            
            vv(ii,kk) = AA(ii,ii)*lam(ii,kk,tt);

            for jj = N_ii
            	vv(ii,kk) = vv(ii,kk) + AA(ii,jj)*lam(jj,kk,tt);
            end
        end
    end
    
    % multiplier average
    lamAvg = mean(lam(:,:,tt));

  % cost updating
    for ii=1:NN
        
        % data extraction
        ci=cc((ii-1)*nni+1 : ii*nni);
        Hi=HH(:,nni*(ii-1)+1:nni*ii);
        UBi=UBB((ii-1)*nni+1:ii*nni);
        LBi=LBB((ii-1)*nni+1:ii*nni);
        Di=DD(ii,SS*(ii-1)+1:SS*ii);
            
        % local solution
        [ZZ(ii,:,tt),~,~]=linprog(ci+(vv(ii,:))*Hi',Di,dd(ii),...
            [],[],LBi,UBi,lpoptions);
        
            % Running average
        if tt==1
            ZRA(ii,:,tt) = ZZ(ii,:,tt);
        else
            ZRA(ii,:,tt) = (1/tt)*((tt-1)*ZRA(ii,:,tt-1)...
                                             +ZZ(ii,:,tt));
        end
        
        ffii = ci * ZZ(ii,:,tt)';
        primCost(tt) = primCost(tt) + ffii;
        %
        ffii = ci * ZRA(ii,:,tt)';
        primRA(tt) = primRA(tt) + ffii;
        %
        qqi = ci * ZZ(ii,:,tt)'...
            +lam(ii,:,tt) * ( ZZ(ii,:,tt) - bi' )';
        dualCost(tt) = dualCost(tt) + qqi;
        %
        qqi = ci *ZRA(ii,:,tt)'...
            + lam(ii,:,tt) * ( ZRA(ii,:,tt) - bi' )';
        dualRA(tt) = dualRA(tt) + qqi;
            
        consErr(:,:,tt) = consErr(:,:,tt) + abs(lam(ii,:,tt) - lamAvg);
    end
    fprintf((repmat('\b', 1,msglnt1+msglnt2)));
    fprintf("\n");
end

%% theoretical basis

% primal problem
% min{z€Z} f(z)
%   == min{z1...zN} sum{1,N}(ci'*zi)
%
% subj: 	sum{1,N} (Hi*zi)=b = 0
%               == h_l(z)=0, l={1...m}
%
%                   Hi*zi-bi==gi(zi)
%
%           ( g_j(z) =< o j€{1...r} )
%
%           zi € Pi, i€{1...N}
%
%### need to implement the Pi containing the limits on Di*zi<=di
%
%-------
%
% lagr L(z,lam)= f(z) + sum{1,m}( lam_l*h_l(z) ) =< f(z)
% 
% grad{z}(L(z*,lam* (,mu*) )=0
%
% ( mu*_j>=0
%   mu*_j*g_j(z*)=0 j={1...r} )
%
%   foreach lam
%   
%------
%
% q(lam (,mu) ) = inf{z€Z}( L )
%
% Dom(q)= {(lam (,mu) ) | q > -oo}
%
%__________________
% dual problem
%
% max{(lam (,mu))€Dom(q) } (q) == min{(lam (,mu))€Dom(q) } (-q) subj:  ( mu
% >=0 )
% 
% q*= sup{lam (,mu>=0) } ( q ) duality: q*=<f* (strong d. if = )
% 
% q* = sup{lam (,mu>=0) }(q) =< inf {z€Z,h(z)=0 (,g(z)=<0)} (f(z)) = f*
%
%_____________________
% distributed algorithm for constraint coupled problems
%
% initialize:    mu_i^0 >=0
%		 any la_i^0
% 
% evolution:	for t>=0
%     gather      mu_j^t from neigh | j€N_i
%       
%     compute     v_i^t+1= sum{j € N_i} (a_ij*mu_j^t)
%     
%                 z_i^t+1 € argmin{z_i€Z_i}(f_i(z_i)+(v_^t+1)'*g_i(z_i))
%                 
%     update      mu_i^t+1 = P{mu>=0} (v_i^t+1+alpha^t*g_i(z_i^t+1))     
%
