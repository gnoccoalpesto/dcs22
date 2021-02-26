function [AA_NW, AA] = pathGraph(NN, stoch)
% pathGraph  Generate a NxN path graph
%   Inputs: N Number of vertexes
%           stoch row/doubly stochasticity
%
%   Outputs: AA_NW Adjacency matrix
%            AA    Weight matrix
%
%   [AA_NW, AA] = pathGraph(p,NN,'row') returns a row stochastic 
%                 weight matrix.
%
%   [AA_NW, AA] = pathGraph(p,NN,'doubly') returns a doubly stochastic 
%                 weight matrix.
%

adj_vect = zeros(NN,1);
adj_vect(2) = 1;
AA_NW = toeplitz(adj_vect);

AA = weightMatrix(AA_NW, stoch);

end