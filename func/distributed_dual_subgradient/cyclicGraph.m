function [AA_NW, AA] = cyclicGraph(NN, stoch)
% cyclicGraph  Generate a NxN cyclic graph
%   Inputs: N Number of vertexes
%           stoch row/doubly stochasticity
%
%   Outputs: AA_NW Adjacency matrix
%            AA    Weight matrix
%
%   [AA_NW, AA] = cyclicGraph(p,NN,'row') returns a row stochastic 
%                 weight matrix.
%
%   [AA_NW, AA] = cyclicGraph(p,NN,'doubly') returns a doubly stochastic 
%                 weight matrix.
%

v = ones(1, NN - 1);
AA_NW = diag(v, 1);
AA_NW(end, 1) = 1;
AA_NW = or(AA_NW,AA_NW'); % Symmetric
% Weight Matrix
AA = weightMatrix(AA_NW, stoch);

end