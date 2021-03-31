% WEIGHTS ci of linear problem
% min{z1...zN} sum{1,N}(ci'*zi)
% subj: 	
% LINEAR COUPLING CONSTRAINT:
%      sum{1,N} (Hi*zi)=b == sum{1,N}(Hi*zi-bi)=0
%      Hi*zi-bi==gi(zi)
% POLYHEDRON DEFINITION Pi
%       zi € Pi, i€{1...N}


% WEIGHTED ADJACENCY MATRIX
%p = 0.1;[AA_NW, AA] = binomialGraph(1, NN, 'doubly');
% by deafult AA: same col/row of 0.2
%            AA_NW missing self loops

%MAXITERS = 5e4;

%ZZ = zeros(NN,MAXITERS);
%ZZ_RA = zeros(NN,MAXITERS);%running average

%LA = zeros(NN,MAXITERS);
%VV = zeros(NN,1);

%primal_cost = zeros(MAXITERS,1);
%dual_cost   = zeros(MAXITERS,1);
%primal_cost_RA = zeros(MAXITERS,1);
%dual_cost_RA   = zeros(MAXITERS,1);

%consensus_err = zeros(MAXITERS,1);

%for tt = 1:MAXITERS-1
   
  %gamma_t = 0.1*(1/tt)^0.6; % 1/tt^alpha with alpha in (0.5, 1]

  %for ii=1:NN
   %N_ii = find(AA_NW(:,ii) == 1)';
   % returns each row element position which is not a self loop, f.e. row

    % update: averaging
    % vi =  sum_{k:1:N} aik*muk^t
    %VV(ii) = AA(ii,ii)*LA(ii,tt);
    %for jj = N_ii % like py: for in [array]
      %VV(ii) = VV(ii) + AA(ii,jj)*LA(jj,tt);
    %end
    % but why "initializing it with the self loop element?
  %end

% Primal Update: local minimization  
% zi^t+1=argmin{zi€Pi}(f(zi) + (vi^t+1)' * gi(zi))
%       =argmin{zi€Pi}( (ci'+ (vi^t+1)'*Hi)*zi)

%   for ii=1:NN
%      XX(ii,tt) = linprog(cc_LP(ii)+VV(ii)*AA_LP(ii),...
%             [],[],[],[],LB(ii),UB(ii),options);
%   end

% ?????????????????????????????????????????????????????????????????????????
%   for ii=1:NN
%     if cc_LP(ii)+VV(ii)*AA_LP(ii)>=0
%         ZZ(ii,tt) = LB(ii);
%     else
%         ZZ(ii,tt) = UB(ii);
%     end
%   end
  
  % Running average
%   for ii=1:NN
%       if tt==1
%           ZZ_RA(ii,tt) = ZZ(ii,tt);
%       else
%           ZZ_RA(ii,tt) = (1/tt)*((tt-1)*ZZ_RA(ii,tt-1)+ZZ(ii,tt));
%       end
%   end
  
  % Dual Update
  
  % mu^t+1=max{0, vi^t+1 + alfa^t *(Hi*zi^t+1 -bi) }
  % mu^{t+1} = mu^t + gamma^t* ( sum_i grad_q_i(mu^t) )

%   for ii=1:NN
%     grad_ii = ZZ(ii,tt)-bb_LP(ii);
% % why like this??????????????????????????????????????????????????????????
%     MU(ii,tt+1) = VV(ii) + gamma_t*grad_ii; % update
%     MU(ii,tt+1) = max(MU(ii,tt+1),0); % projection   
%   end
  
  % Performance check ?????????????????????????????????????????????????????
% 	MU_avg = mean(MU(:,tt));

%   for ii=1:NN
%     ff_ii = cc_LP(ii)*ZZ(ii,tt);
%     primal_cost(tt) = primal_cost(tt) + ff_ii;
%     
%     ff_ii = cc_LP(ii)*ZZ_RA(ii,tt);
%     primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
%     
%     qq_ii = cc_LP(ii)*ZZ(ii,tt) + MU(ii,tt)*(ZZ(ii,tt)-bb_LP(ii));
%     dual_cost(tt) = dual_cost(tt) + qq_ii;
%     
%     qq_ii = cc_LP(ii)*ZZ_RA(ii,tt) + MU(ii,tt)*(ZZ_RA(ii,tt)-bb_LP(ii));
%     dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;
% % what's the difference between RA and not???????????????????????????????
%
%     consensus_err(tt) = consensus_err(tt) + norm(MU(ii,tt) - MU_avg);
%   % why uses norm (distance actual-average)??????????????????????????????
%   end
% end
%%
% Last value [for plot]
% tt = MAXITERS;
% fprintf('Iteration n. %d\n',tt);
% MU_avg = mean(MU(:,tt));
% average @ last+1 iterations

% for ii=1:NN
%   N_ii = find(AA_NW(:,ii) == 1)';
%   % extract rows element not zero(i.e. the self loops)
%   
%   % v_i =  sum_{j in N_i U { i } } aw_ij mu_j^t
%   % aw_ij elements of weig.adj.mat.
%   VV(ii) = AA(ii,ii)*MU(ii,tt);
%   
%   for jj = N_ii% py for x in y
%     VV(ii) = VV(ii) + AA(ii,jj)*MU(jj,tt);
%   end
% %since ii,ii always zero, why separate the calcolus? it should have sense
% % if first part is calculated on find(AA_NW(:,ii) != 1)'?????????????????
% end

% for ii=1:NN
%     % what's this branching??????????????????????????????????????????????
%     if cc_LP(ii)+VV(ii)*AA_LP(ii)>=0
%         ZZ(ii,tt) = LB(ii);
%     else
%         ZZ(ii,tt) = UB(ii);
%         % and why selects only @ bounds ?????????????????????????????????
%     end
%     ZZ_RA(ii,tt) = (1/tt)*((tt-1)*ZZ_RA(ii,tt-1)+ZZ(ii,tt));
% 
%     ff_ii = cc_LP(ii)*ZZ(ii,tt);
%     primal_cost(tt) = primal_cost(tt) + ff_ii;
%     
%     ff_ii = cc_LP(ii)*ZZ_RA(ii,tt);
%     primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
%     
%     qq_ii = cc_LP(ii)*ZZ(ii,tt) + MU(ii,tt)*(ZZ(ii,tt)-bb_LP(ii));
%     dual_cost(tt) = dual_cost(tt) + qq_ii;
%     
%     qq_ii = cc_LP(ii)*ZZ_RA(ii,tt) + MU(ii,tt)*(ZZ_RA(ii,tt)-bb_LP(ii));
%     dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;
%     % as before
% end

%%
figure
% why the cost is negative in first instance???????????????????????????????
  semilogy(1:MAXITERS,abs(primal_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogy(1:MAXITERS,abs(primal_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
  semilogy(1:MAXITERS,abs(dual_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  semilogy(1:MAXITERS,abs(dual_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('cost error')
  legend('primal cost','primal cost with running avg','dual cost',...
      'dual cost with running avg')
  
%%
figure
  plot(1:MAXITERS,MU, 'LineWidth',2);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('\mu_i^t')

%%
figure
  plot(1:MAXITERS,sum(ZZ,1)-bb_centr, 'LineWidth',2);
  hold on, grid on, zoom on
  plot(1:MAXITERS,sum(ZZ_RA,1)-bb_centr, 'LineWidth',2);
  xlabel('t')
  ylabel('x_1^t + ... + x_N^t - b')
  legend('x','x from running avg')
  
  
%%
figure
  semilogy(1:MAXITERS,consensus_err(1:MAXITERS), 'LineWidth',2);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('consensus error')
  
  %% conclusions
  % dual cost: converging rapidly (orange), then oscillate
  % (violet) update on running average
  % (blue plot) zi^t non converging, oscillating
  % (red) running average is converging
  %
  % convergence on the mu variables
  % convergence of running averages
