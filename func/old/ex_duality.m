% exemple on duality and primal problem

% primal:
% min{z} (c' * z)
% subj: a_l' * z=b_l l={1...m}
%       z>=0
%       
% L= c' * z +sum{1,m}(lam_l(b_l-a_l' * z)
%  = (c-sum{1,m}(lam_l * a_l) )' * z + sum{1,m}(lam_l * b_l)
%  
% foreach z>=0 : inf @ z=0 if c-sum{1,m}(lam_l * a_l) >=0
% q(lam)= sum{1,m}(lam_l * b_l)
% Dom(q)={lam € R^m | c-sum{1,m}(lam_l * a_l) >=0 }
%  
% dual:
% max{lam}(sum{1,m}(lam_l * b_l))
% subj: sum{1,m}(lam_l * a_l_j) =< c_j, j={1...d}
% 
% subgrad alg:
% only subj: mu>=0 ??????
% 
% lam^t+1= lam^t + alfa^t * subgrad{lam}(q (lam^t)) , no projections
% (for mu^t+1: same, with proj{mu>=0}== max(mu^t+1,0))
% {alfa^t}{t>=0} stepsize (usual rules)
% 
% subgrad{lam}(q)=[...h_l(z^t+1))...]'
% 
% z^t+1=argmin{z€Z}(L (z, lam^t) )
% lam_t+1= lam^t + alfa^t * h(z^t+1)
