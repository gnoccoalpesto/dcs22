
	%clear all; close all; clc;
	N = 3; % NUMBER OF AGENTS 
	M = 4; % SIZE OF AGENTS' LOCAL OPTIM. VAR. ATTENTION! IT MUST STAND M>N

	%%% MODEL A
	% % profit associated with assigning job j to machine i
	P_ = randi([10 25],N,M);
	P  = max(P_,[],'all') +1 - P_ ;

	% claim on the capacity of machine i by job j
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

	H = repmat(eye(M,M), [1,N]); % == [...Hj...]
	b = ones(M,1); 
	% nota: la "separazione" della b sui vari agenti è arbitraria:
	% - tutta ad uno
	% - media a ciascuno
	% - whatever ...

	% machine capacity, local constraints
	D = zeros(N,M*N);% D=[...Dj...]
	for i = 1:N
	    D(i,M*(i-1)+1:M*(i-1)+M) = W(i,:);
	end

	%maximize profit
	c = reshape(P',[1,M*N]);% c=[...cj...]
	LB = zeros((N*M),1);
	UB = ones((N*M),1);

	%[zopt, fopt] = linprog(c,D,d,H,b,LB,UB)
	% min{x} c'*z | A*z=<b & Aeq*z = beq & lb=<z=<ub

	clearvars -except c D d H b LB UB zopt fopt