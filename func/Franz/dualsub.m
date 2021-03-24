function...
    [primCost,dualCost,primRA,dualRA,...consErr,
                lam,ZZ,ZRA,MAXITERS]=...
                           dualsub(NN,nni,AA,AANW,bb,cc,dd,DD,HH,LBB,UBB)
% dd,DD,LBB unused
    MAXITERS = 1e5 ;
    fprintf('iterations: %d% \n',MAXITERS);

    fprintf('initializing...');
    primCost=zeros(MAXITERS,1); % cost € R
    dualCost=zeros(MAXITERS,1);
    primRA=zeros(MAXITERS,1);
    dualRA=zeros(MAXITERS,1);
%     consErr=zeros(MAXITERS,1);
    
    % agents init
    % agent zi € R^ni, zz=[...;zi;...] agent matrix € R^N x ni
    % ZZ=[zz_1st_iteraction, ... , zz_MAXITERS]
    ZZ=zeros(NN,nni,MAXITERS); 
    ZRA=zeros(NN,nni,MAXITERS);
    
    SS=nni;
    % when %==SS is displayed, means that it's our particular case nni==SS
    
    %multiplier init
    lam=zeros(NN,SS,MAXITERS);% ????????????????????????????????????????
    vv=zeros(NN,SS); % vv=[...;vi;...], vi€ R ^ SS
    
    % b splitting
    bi=bb/NN;
    
    for tt = 1:MAXITERS-1
        % todo
        %if mod(tt,500)==0,fprintf('progress: %d\n',tt*500/MAXITERS);end
        
        %stepsize
        alpha_t = 0.1*(1/tt)^.6; % 1/tt^alpha with alpha in (0.5, 1]
        
        
        
        %gather la_j^t from neigh | j€N_i
        %moved after lamda init
        
        %compute	vm_i^t+1= sum{j € N_i} (a_ij*mu_j^t)
	%		vl_i^t+1= sum{j € N_i} (a_ij*la_j^t)
            % == update: averaging
            % vi =  sum_{k:1:N} aik*muk^t
            
%         for ii=1:NN %modified neigh version
%             %-2-
%             N_ii = find(AA(:,ii)~=0)'; % use eps not 0
%             for jj = N_ii,vv(ii) = vv(ii) + AA(ii,jj)*lam(jj,tt);end
%         end

%         for ii=1:NN % full row version, no for kk=1:SS
%             N_ii = find(AANW(:,ii) == 1)';
%             vv(ii,:) = AA(ii,ii)*lam(ii,:,tt);
%             for jj = N_ii,vv(ii) = vv(ii) + AA(ii,jj)*lam(jj,tt);end
%         end

        for ii=1:NN %iterate for each agent
            for kk=1:SS %iterate in agent's dimension
                N_ii = find(AANW(:,ii) == 1)';% neigh(i)

        % update: averaging
        % vvi =  sum_{k:1:N} aik*lam_k^t
                vv(ii,kk) = AA(ii,ii)*lam(ii,kk,tt);

                for jj = N_ii
                    vv(ii,kk) = vv(ii,kk) + AA(ii,jj)*lam(jj,kk,tt);
                end
            end
        end
        
        
    
        %        z_i^t+1 € argmin{z_i€Z_i}(f_i(z_i)+(v_^t+1)'*g_i(z_i))
        %             zi^t+1=argmin{zi€Pi}(ci*zi + (vi^t+1)' * gi(zi))
        %            =argmin{zi€Pi}( (ci'+ (vi^t+1)'*Hi)*zi)
        

        
        for ii=1:NN
            
            ci=cc((ii-1)*nni+1 : ii*nni);
            Hi=HH(:,nni*(ii-1)+1:nni*ii);
            UBi=UBB((ii-1)*nni+1:ii*nni);
            
            %%%         for ii=1:NN, H(:,nni*(ii-1)+1:nni*(ii-1)+nni);end     %%%
%                           ==Hii € R^ SS x nni
%         to select only the matrix Hi, relatively to agent ii
%           takes all rows and only the columns spanning for 
%            agent's dimension
            ZZ(ii,:,tt)=linprog(ci',[],[],Hi,vv(ii,:),[],UBB);
            % introdurre vincolo f(xi)+(v^t+1)'*g(xi)
            % Running average
            for jj=1:nni %==SS
               % -a- 
                if tt==1
                    ZRA(ii,jj,tt) = ZZ(ii,jj,tt);
                else
                    ZRA(ii,jj,tt) = (1/tt)*((tt-1)*ZRA(ii,jj,tt-1)...
                                        +ZZ(ii,jj,tt));
                end
            end
        end
    
        
        %update      mu_i^t+1 = P{mu>=0} (vm_i^t+1+alpha^t*g_i(z_i^t+1))
	%		la_i^t+1=v_i^t+1+alpha^t*g_i(z_i^t+1)
    
        for ii=1:NN
            
            Hi=HH(:,nni*(ii-1)+1:nni*ii);
            for jj=1:SS
                grad_ii_jj= Hi(jj,jj)* ZZ(ii,jj,tt)-bi(jj);% -4-
                        
            lam(ii,jj,tt+1) = vv(ii,jj) + alpha_t*grad_ii_jj;
            end
        end
               

         % Running average
%         for ii=1:NN
%             if tt==1
%                 zra(ii,tt) = z(ii,tt);
%             else
%                 zra(ii,tt) = (1/tt)*((tt-1)*zra(ii,tt-1)+z(ii,tt));
%             end
%         end
        
    
%         lamAvg=mean(lam(:,:,tt)); %paused

        for ii=1:NN
            ci=cc((ii-1)*nni+1 : ii*nni);
            
            for jj=1:nni % == SS
                
            ffii_jj = ci(jj) * ZZ(ii,jj,tt);
            
            primCost(tt) = primCost(tt) + ffii_jj;
            %
            ffii_jj = ci(jj) * ZRA(ii,jj,tt);
            
            primRA(tt) = primRA(tt) + ffii_jj;
            %
            qqi = ci(jj) * ZZ(ii,jj,tt)...
                ...
                +lam(ii,jj,tt) * ( ZZ(ii,jj,tt) - bi(jj) );
            
            dualCost(tt) = dualCost(tt) + qqi;
            %
            qqi = ci(jj) *ZRA(ii,jj,tt)...
                ...
                + lam(ii,jj,tt) * ( ZRA(ii,jj,tt) - bi(jj) );
            
            dualRA(tt) = dualRA(tt) + qqi;

%             consErr(tt) = consErr(tt) + norm(lam(ii,tt) - lamAvg);
            end
        end
    
    end
    
    tt=MAXITERS;
%     lamAvg=(1/tt)*((tt-1)*lamAvg+lam(ii,tt)); %"running" avg
    
    for ii=1:NN %iterate for each agent
        for kk=1:SS %iterate in agent's dimension
            N_ii = find(AANW(:,ii) == 1)';% neigh(i)

        % update: averaging
        % vvi =  sum_{k:1:N} aik*lam_k^t
            vv(ii,kk) = AA(ii,ii)*lam(ii,kk,tt);

            for jj = N_ii
            	vv(ii,kk) = vv(ii,kk) + AA(ii,jj)*lam(jj,kk,tt);
            end
        end
    end
    
  
    for ii=1:NN
        
        ci=cc((ii-1)*nni+1 : ii*nni);
        Hi=HH(:,nni*(ii-1)+1:nni*ii);
        UBi=UBB((ii-1)*nni+1:ii*nni);
        
        for jj=1:nni % == SS
        
        if ci(jj)+ Hi(jj,jj) * vv(ii,jj)<0
               
        	ZZ(ii,jj,tt) = UBi(jj);
        end
        
        % Running average
        ZRA(ii,jj,tt) = (1/tt)*((tt-1)*ZRA(ii,jj,tt-1)...
                                        +ZZ(ii,jj,tt));
        
        
        

        ffii_jj = ci(jj) * ZZ(ii,jj,tt);

        primCost(tt) = primCost(tt) + ffii_jj;
        %
        ffii_jj = ci(jj) * ZRA(ii,jj,tt);

        primRA(tt) = primRA(tt) + ffii_jj;
        %
        qqi = ci(jj) * ZZ(ii,jj,tt)...
            ...
            +lam(ii,jj,tt) * ( ZZ(ii,jj,tt) - bi(jj) );

        dualCost(tt) = dualCost(tt) + qqi;
        %
        qqi = ci(jj) *ZRA(ii,jj,tt)...
            ...
            + lam(ii,jj,tt) * ( ZRA(ii,jj,tt) - bi(jj) );

        dualRA(tt) = dualRA(tt) + qqi;

%             consErr(tt) = consErr(tt) + norm(lam(ii,tt) - lamAvg);
        end
    end
    
     
end

%% comments
% xi^t (blue one) oscillates, not converging
% running avarage converges 
%       XX_RA(ii,tt) = (1/tt)*((tt-1)*XX_RA(ii,tt-1)+XX(ii,tt))
% lam^t converging as expected for the dual update
%
%% primal problem
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