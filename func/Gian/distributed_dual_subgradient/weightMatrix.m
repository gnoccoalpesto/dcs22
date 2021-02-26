function AA = weightMatrix(AA_NW, stoch)
% weightMatrix  Generate a NxN weight matrix associated to an adjacency
% matrix
%   Inputs: AA_NW Adjacency matrix
%           stoch row/doubly stochasticity
%
%   Outputs: AA Weight matrix
%
%   [AA_NW, AA] = weightMatrix(p,NN,'row') returns a row stochastic 
%                 weight matrix.
%
%   [AA_NW, AA] = weightMatrix(p,NN,'doubly') returns a doubly stochastic 
%                 weight matrix.
%
NN = size(AA_NW,1);
AA = zeros(NN,NN);

DEGREE = sum(AA_NW);
if strcmp(stoch,'row')
    for ii = 1:NN
      N_ii = find(AA_NW(ii,:) == 1)';
      for jj = [N_ii; ii]
		AA(ii,jj) = 1/(1 + DEGREE(ii)); % row stoch
      end
    end
elseif strcmp(stoch,'doubly')
    for ii = 1:NN
          N_ii = find(AA_NW(:,ii) == 1)';
       for jj = N_ii
         AA(ii,jj) = 1/(1 + max(DEGREE(ii),DEGREE(jj) )); % doubly stoch
       end
       AA(ii,ii) = 1 - sum(AA(ii,:));
    end
else
	error('stoch must be a string equal to row or doubly')
end

end