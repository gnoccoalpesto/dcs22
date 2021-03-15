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


%%Primal problem

MAXITERS = 5e4;

ZZ = zeros(NN,MAXITERS); 
ZZ_RA = zeros(NN,MAXITERS);

LA = zeros(MAXITERS,1);
%H=eye(M);

primal_cost = zeros(MAXITERS,1);
dual_cost   = zeros(MAXITERS,1);
primal_cost_RA = zeros(MAXITERS,1);
dual_cost_RA   = zeros(MAXITERS,1);

%che è sta roba?
% tt_ra = 1;
% NMIN=1;

for tt = 1:MAXITERS-1
    1==1;
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