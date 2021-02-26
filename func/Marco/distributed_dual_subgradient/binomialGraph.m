function [AA_NW, AA] = binomialGraph(p,NN,stoch)
% binomialGraph  Generate a NxN random graph with edge probability p
%   Inputs: p Edge probability
%           N Number of vertexes
%           stoch row/doubly stochasticity
%
%   Outputs: AA_NW Adjacency matrix
%            AA    Weight matrix
%
%   [AA_NW, AA] = binomialGraph(p,NN,'row') returns a row stochastic 
%                 weight matrix.
%
%   [AA_NW, AA] = binomialGraph(p,NN,'doubly') returns a doubly stochastic 
%                 weight matrix.
%


I_NN = eye(NN,NN);
notI = ~I_NN;

while 1
  AA_NW = binornd(1,p,NN,NN); % Randomly generated matrix
  % Adj = rand(1,p,NN,NN);

  AA_NW = AA_NW.*notI; % Set diagonal to 0 (remove self-Loops)
  AA_NW = or(AA_NW,AA_NW'); % Symmetric

  % Test connection
  test = (I_NN+AA_NW)^NN;
  if ~any(any(~test))
    break
  end  
end

% Weight Matrix
AA = weightMatrix(AA_NW, stoch);
end