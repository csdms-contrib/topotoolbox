function c = getcovariate(P,cova)

%GETCOVARIATE reads covariate nal from P
%
% Syntax
%
%     c = getcovariate(P,cova);
%
% Description
%
%     GETCOVARIATE tries to make sense of a covariate supplied by the user.
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019
    

% get covariate
if isempty(cova)
    
    % With no covariate, the function will return the distance upstream
    c = P.S.distance;
    
elseif iscell(cova)
    
    % cell array
    c = cellfun(@(cc) getcovariate(P,cc),cova,'UniformOutput',false);
    c = horzcat(c{:});
    return

elseif ischar(cova)
            
    if strcmp(cova,'z')
        c = P.z;
    elseif strcmp(cova,'x')
        c = P.S.x;
    elseif strcmp(cova,'y')
        c = P.S.y;
    elseif strcmp(cova,'distance')
        c = P.S.distance;
    else
        try 
            c = P.covariates.(cova);
        catch
            error('Cannot interpret string');
        end
    end

elseif isa(cova,'GRIDobj')
    
    c = getnal(P.S,cova);
    
elseif size(cova,1) == numel(P.S.x)
    
    % cova is one or several node attribute lists
    c = cova;
    
elseif istable(cova)
    
    if size(cova,1) ~= numel(P.S.x)
        error(['Table must have ' num2str(numel(P.S.x)) ' rows.']);
    end
    
    
else
    
    error('Cannot handle input.');
    
end