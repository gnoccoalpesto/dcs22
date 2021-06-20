% spawns entities in a R^n (R^2) space
% confinement area can be a circle or a rectangle
% input: 
%        togglePlot (bool)
%        nttNum: number of entities to spawn
%        shape:{square, rectangle} determined by number of given dim
%
% output: returnVector: coordinates of spawned entities
%
function [returnVector] =spawnEntities (togglePlot,nttNum,spnAreaWidth,...
                                    spnAreaHeight,rngSeed)
    if nargin>4
        if rngSeed==0, rng('shuffle');% time seeding
        else, rng(rngSeed);
        end
    else, rng(1); 
    end
%NOTE usage: use the same value for rng to replicate results

    envBaseDim=2; % environment base dimension
%NOTE practilcally unused, may generate in 3 dimensions too

% result init
    returnVector=zeros(nttNum,envBaseDim);
    % result randomization
    randVector=rand(size(returnVector));
    
    %result scaling
    if nargin<3% unitary square
        returnVector=randVector;
    elseif nargin<4% circular shape
%         
        returnVector=randVector.*[spnAreaWidth/2,2*pi];
        returnVector=[returnVector(:,1).*cos(returnVector(:,2)),...
             returnVector(:,1).*sin(returnVector(:,2))]+spnAreaWidth/2;
%          remove +spawnAreaDim to center spawn area center as origin
    else% rectangle
        returnVector=randVector.*[spnAreaWidth spnAreaHeight];
%        add -[spnAreaWidth/2 spnAreaHeight/2] to center on origin
    end
    
% display the plot    
    if togglePlot
        plot(returnVector(:,1),returnVector(:,2),'*')
    end
end