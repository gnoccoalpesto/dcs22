function [AA_NW, AA] = myGraph(NN, stoch, tyGraph,prob)
% myGraph calls the adequate function to generate agents graph
%   Inputs: N Number of verteces
%           'row','doubly' stochasticity
%	    type of generated graph: 
%		'cyclic'
%		'path'
%		'binomial'	
%
%   Outputs: AA_NW Adjacency matrix
%            AA    Weight matrix
%
%NOTE function are provided by A.Testa


    if strcmp(tyGraph,'cyclic')
        [AA_NW,AA]=cyclicGraph(NN,stoch);
    elseif strcmp(tyGraph,'path')
        [AA_NW,AA]=pathGraph(NN,stoch);
    elseif strcmp(tyGraph,'binomial') && nargin>3 && abs(prob)<=1
        [AA_NW,AA]=binomialGraph(abs(prob),NN,stoch);
    else
        [AA_NW,AA]=binomialGraph(1,NN,stoch);
    end

end
