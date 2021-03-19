function...
    [primCost,dualCost,primRA,dualRA,consErr,lam,zz,zra,MAXITERS]=...
                           dualsub(NN,AA,AANW,bb,cc,dd,DD,HH,LBB,UBB)
% dd,DD,LBB unused
    MAXITERS = 5e4;
    fprintf('progress %d%\n',MAXITERS);

    disp('initializing...');
    primCost=zeros(MAXITERS,1);
    dualCost=zeros(MAXITERS,1);
    primRA=zeros(MAXITERS,1);
    dualRA=zeros(MAXITERS,1);
    consErr=zeros(MAXITERS,1);
    
    zz=zeros(NN,MAXITERS); 
    zra=zeros(NN,MAXITERS);
    
    %multiplier init
    lam=zeros(MAXITERS,1);
    vv=zeros(NN,1);
    
    for tt = 1:MAXITERS-1
        if mod(tt,500)==0,fprintf('progress: %d\n',tt*.02);end
        
        %stepsize
        alpha_t = 0.1*(1/tt)^.6; % 1/tt^alpha with alpha in (0.5, 1]
        
        
        
        %gather la_j^t from neigh | j€N_i
        %moved after lamda init
        
        %compute	vm_i^t+1= sum{j € N_i} (a_ij*mu_j^t)
	%		vl_i^t+1= sum{j € N_i} (a_ij*la_j^t)
            % == update: averaging
            % vi =  sum_{k:1:N} aik*muk^t
            
%         for ii=1:NN
%             %-2-
%             N_ii = find(AA(:,ii)~=0)';
%             for jj = N_ii
%                 vv(ii) = vv(ii) + AA(ii,jj)*lam(jj,tt);
%             end
%         end
        for ii=1:NN
            N_ii = find(AANW(:,ii) == 1)';

    % update: averaging
    % vi =  sum_{k:1:N} aik*muk^t
            vv(ii) = AA(ii,ii)*lam(ii,tt);
    
            for jj = N_ii
                vv(ii) = vv(ii) + AA(ii,jj)*lam(jj,tt);
            end
        end
        
    
        %        z_i^t+1 € argmin{z_i€Z_i}(f_i(z_i)+(v_^t+1)'*g_i(z_i))
        %             zi^t+1=argmin{zi€Pi}(ci*zi + (vi^t+1)' * gi(zi))
        %            =argmin{zi€Pi}( (ci'+ (vi^t+1)'*Hi)*zi)
        
        for ii=1:NN
            %-3-
            if cc(ii)+vv(ii)*HH(ii)<0
                zz(ii,tt) = UBB(ii);
            end
            % Running average
            if tt==1
                zra(ii,tt) = zz(ii,tt);
            else
                zra(ii,tt) = (1/tt)*((tt-1)*zra(ii,tt-1)+zz(ii,tt));
            end
        end
        
        %update      mu_i^t+1 = P{mu>=0} (vm_i^t+1+alpha^t*g_i(z_i^t+1))
	%		la_i^t+1=v_i^t+1+alpha^t*g_i(z_i^t+1)
    
         for ii=1:NN
            grad_ii = HH(ii)*zz(ii,tt)-bb(ii);
            lam(ii,tt+1) = vv(ii) + alpha_t*grad_ii;   
         end   
        
         % Running average
%         for ii=1:NN
%             if tt==1
%                 zra(ii,tt) = z(ii,tt);
%             else
%                 zra(ii,tt) = (1/tt)*((tt-1)*zra(ii,tt-1)+z(ii,tt));
%             end
%         end
        
    
        lamAvg=mean(lam(:,tt));
        for ii=1:NN
            ffi = cc(ii)*zz(ii,tt);
            primCost(tt) = primCost(tt) + ffi;

            ffi = cc(ii)*zra(ii,tt);
            primRA(tt) = primRA(tt) + ffi;

            qqi = cc(ii)*zz(ii,tt) + lam(ii,tt)*(zz(ii,tt)-bb(ii));
            dualCost(tt) = dualCost(tt) + qqi;

            qqi = cc(ii)*zra(ii,tt) + lam(ii,tt)*(zra(ii,tt)-bb(ii));
            dualRA(tt) = dualRA(tt) + qqi;

            consErr(tt) = consErr(tt) + norm(lam(ii,tt) - lamAvg);
        end
    
    end
    
    tt=MAXITERS;
    lamAvg=(1/tt)*((tt-1)*lamAvg+lam(ii,tt)); %"running" avg
    
    for ii=1:NN
        N_ii = find(AA(:,ii)~=0)';
        for jj = N_ii
            vv(ii) = vv(ii) + AA(ii,jj)*lam(jj,tt);
        end
    end
    
    for ii=1:NN
        if cc(ii)+vv(ii)*HH(ii)<0
            zz(ii,tt) = UBB(ii);
        end
        % Running average
        zra(ii,tt) = (1/tt)*((tt-1)*zra(ii,tt-1)+zz(ii,tt));
        
        ffi = cc(ii)*zz(ii,tt);
        primCost(tt) = primCost(tt) + ffi;

        ffi = cc(ii)*zra(ii,tt);
        primRA(tt) = primRA(tt) + ffi;

        qqi = cc(ii)*zz(ii,tt) + lam(ii,tt)*(zz(ii,tt)-bb(ii));
        dualCost(tt) = dualCost(tt) + qqi;

        qqi = cc(ii)*zra(ii,tt) + lam(ii,tt)*(zra(ii,tt)-bb(ii));
        dualRA(tt) = dualRA(tt) + qqi;

        consErr(tt) = consErr(tt) + norm(lam(ii,tt) - lamAvg);
        
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