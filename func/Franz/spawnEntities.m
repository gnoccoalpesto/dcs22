% spwns entities in a R^n (R^2) space
% param: number2spawn: number of entities
%        dimensions: vector of space dimensions
%        shapeType: can be square or circular, may add any shape type
% returns: returnVector: coordinates of spawned entities
function [returnVector] =spawnEntities (number2spawn,spawnAreaDim,shapeType)
    % rng init
    rng('shuffle') % uses time to seed
    
    % specify in which R^n entities are spawned
    spaceBaseDim=2;
    % coord vector init
    returnVector=zeros(number2spawn,spaceBaseDim);
    randVector=rand(size(returnVector));
    % checks if shape is specified
    if nargin<3
        % defaut shape: square
        returnVector=randVector.*spawnAreaDim;%-spawnAreaDim/2;
        % remove comment of previous line to move rectangles center to origin
    elseif shapeType=="circular"
        returnVector=randVector.*[spawnAreaDim/2,2*pi];
        returnVector=[returnVector(:,1).*cos(returnVector(:,2)),...
             returnVector(:,1).*sin(returnVector(:,2))]+spawnAreaDim/2;
         % remove +spawnAreaDim to center spawn area center as origin
    end
    togglePlot=1;
    if togglePlot
        plot(returnVector(:,1),returnVector(:,2),'*')
end
