% spawns entities in a R^n (R^2) space
% confinement area can be a circle or a rectangle
% param: 
%        togglePlot (bool)
%        nttNum: number of entities to spawn
%        shape:{square, rectangle} determined by number of given dim
%
% returns: returnVector: coordinates of spawned entities
%
function [returnVector] =spawnEntities (togglePlot,nttNum,spnAreaWidth,...
                                    spnAreaHeight,rngSeed)
    if nargin>4
        if rngSeed==0, rng('shuffle');
        else, rng(rngSeed);
        end
    else, rng(1); 
    end% prefer a known number for replicate a nice configuration
%     
    envBaseDim=2; % environment base dimension
%     
    returnVector=zeros(nttNum,envBaseDim);
    randVector=rand(size(returnVector));
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
%         add: move the space such that its center its origin
    end
%     
    if togglePlot
        plot(returnVector(:,1),returnVector(:,2),'*')
    end
end