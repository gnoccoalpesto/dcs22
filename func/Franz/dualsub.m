    %%
% % Centralized Solution
% 
% % centralized (parallel) dual subgradient:
% % zi^t+1= argmin{zi € Pi}(ci'*zi +(mu^t)'*(Hi*zi-bi)
% %       = argmin{zi € Pi}( ci'+(mu^t)'*Hi)*zi -(mu^t)'*bi)
% % mu^t+1 = max {0, mu^t+alfa^t *sum{1,N}(Hi*zi^t+1 -bi) }
% 
% 
% options = optimoptions('linprog','Display','none');
% [~, fopt, exit_flag] = linprog(c,H,bb_centr,[],[],LB,UB,options);
% 
% if exit_flag ~= 1
%   fprintf(2,'A problem occurred in the centralized solution\n');
%   return;
% end
% 
% fprintf('Centralized optimal cost is %.4g\n',fopt);

% primal problem
% min{z€Z} f(z)
% subj: h_l(z)=0, l={1...m}
%       ( g_j(z) =< o j€{1...r} )
% 
% lagr L(z,lam)= f(z) + sum{1,m}( lam_l*h_l(z) ) =< f(z)
% 
% grad{z}(L(z*,lam* (,mu*) )=0
% ( mu*_j>=0
%   mu*_j*g_j(z*)=0 j={1...r} )
%   foreach lam
%   
% q(lam (,mu) ) = inf{z€Z}( L )
% Dom(q)= {(lam (,mu) ) | q > -oo}
% 
% dual problem
% max{(lam (,mu))€Dom(q) } (q) == min{(lam (,mu))€Dom(q) } (-q)
% subj:  ( mu >=0 )
% 
% q*= sup{lam (,mu>=0) } ( q )
% duality: q*=<f* (strong d. if = )
% 
% q* = sup{lam (,mu>=0) }(q) =< inf {z€Z,h(z)=0 (,g(z)=<0)} (f(z)) = f*


% distributed algorithm for constraint coupled problems

% initialize:    mu_i^0 >=0
% 
% evolution:	for t>=0
%     gather      mu_j^t from neigh | j€N_i
%     
%     compute     v_i^t+1= sum{j € N_i} (a_ij*mu_j^t)
%     
%                 z_i^t+1 € argmin{z_i€Z_i}(f_i(z_i)+(v_^t+1)'*g_i(z_i))
%                 
%     update      mu_i^t+1 = P{mu>=0} (v_i^t+1+gamma^t*g_i(z_i^t+1))
    
    
    
function [primCost,dualCost,primRA,dualRA,consErr]=...
                           dualsub(NN,AA,AANW,bb,cc,dd,DD,HH,LBB,UBB)
    MAXITERS = 5e4;

    primCost=zeros(MAXITERS,1);
    dualCost=zeros(MAXITERS,1);
    primRA=zeros(MAXITERS,1);
    dualRA=zeros(MAXITERS,1);
    consErr=zeros(MAXITERS,1);
    
    ZZ=zeros(NN,MAXITERS); 
    ZZRA=zeros(NN,MAXITERS);
    
    %multiplier init
    LA=zeros(MAXITERS,1);
    VV=zeros(NN,1);
    
    %che è sta roba?
    % tt_ra = 1;
    % NMIN=1;

    for tt = 1:MAXITERS-1
        %gather la_j^t from neigh | j€N_i
        
        %compute	v_i^t+1= sum{j € N_i} (a_ij*mu_j^t)
        %        z_i^t+1 € argmin{z_i€Z_i}(f_i(z_i)+(v_^t+1)'*g_i(z_i))
       
        %update      mu_i^t+1 = P{mu>=0} (v_i^t+1+gamma^t*g_i(z_i^t+1))
        %stepsize
        gamma_t = 0.1*(1/tt)^0.6; % 1/tt^alpha with alpha in (0.5, 1]
        
    end
end

%% comments
%this is to implement an inequality constraint
% MU(tt+1) = MU(tt);
%   for ii=1:NN
%     grad_ii = ZZ(ii,tt)-b(ii);
%     
%     MU(tt+1) = MU(tt+1) + alpha_t*grad_ii;
%   end
% 	MU(tt+1) = max(0,MU(tt+1)); %projection

% xi^t (blue one) oscillates, not converging
% running avarage converges 
%       XX_RA(ii,tt) = (1/tt)*((tt-1)*XX_RA(ii,tt-1)+XX(ii,tt))
% mu^t converging as expected for the dual update