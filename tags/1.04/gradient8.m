function [G,ASP] = gradient8(DEM,res,unit)

% 8-connected neighborhood gradient and aspect of a digital elevation model
%
% Syntax
%
%     [G,ASP] = gradient8(DEM,cellsize,unit)
%
% Description
%
%     gradient8 returns the numerical steepest downward gradient and aspect 
%     of a digital elevation model using an 8-connected neighborhood. 
%
% Input
%
%     DEM       m x n matrix with elevation values
%     cellsize  horizontal spacing between data points (default = 1)
%     unit      'tan' --> tangent (default)
%               'rad' --> radians
%               'deg' --> degrees
%               'sin' --> sine
% 
% Output
%
%     G         matrix with tangent of the gradient
%     ASP       matrix with aspect containing integers from 1 to 8
%               according to the direction of the slope counted clockwise
%               from top.
% 
%                8 1 2
%                7 c 3
%                6 5 4
% 
%               ASP is zero for cells without downward neighbor. 
%                  
% Example
% 
%     [X,Y,DEM] = peaks(100);
%     res = abs(X(1,1)-X(1,2));
%     [G,ASP] = gradient8(DEM,res);
%     subplot(1,2,1)
%     surf(X,Y,DEM,G); shading interp; camlight
%     subplot(1,2,2)
%     surf(X,Y,DEM,ASP); shading flat; camlight
%
%
% See also: EZFLOWACC, CURVATURE, ASPECT
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 22. July, 2010


if nargin == 1;
    res = 1;
    unit = 'tan';
elseif nargin == 2 || nargin == 3;
    if isempty(res)
        res = 1;
    else
        validateattributes(res, {'numeric'}, {'>',0,'scalar'},'gradient8', 'res')
    end
    if nargin == 2
        unit = 'tan';
    end
end


if nargout == 1;
    % check for nans;
    I = isnan(DEM);
    flagnan = any(I(:));
    if flagnan
        DEM(I) = inf;
    end
    
    NEIGH = false(3);
    
    % calculate along orthogonal neighbors
    NEIGH(2,:) = true;
    NEIGH(:,2) = true;
    G = (DEM-imerode(DEM,NEIGH))/res;
    
    % calculate along diagonal neighbors
    NEIGH(:,:) = false;
    NEIGH([1 5 9 3 7]) = true;
    G = max(G,(DEM-imerode(DEM,NEIGH))/hypot(res,res));
    
    if flagnan
        G(I) = nan;
    end
    
elseif nargout == 2;
    
    siz    = size(DEM);
    % pad dem with nans
    DEM    = [nan(1,siz(2)+2); ...
        [nan(siz(1),1) DEM nan(siz(1),1)]; ...
        nan(1,siz(2)+2)];
    
    % function handles for shifting DEM
    neighfun = cell(8,2);
    res2 = hypot(res,res);
    
    neighfun{1,1} = @(x) (x(2:end-1,2:end-1)-x(1:end-2,2:end-1))/res;
    neighfun{2,1} = @(x) (x(2:end-1,2:end-1)-x(1:end-2,3:end))/res2;
    neighfun{3,1} = @(x) (x(2:end-1,2:end-1)-x(2:end-1,3:end))/res;
    neighfun{4,1} = @(x) (x(2:end-1,2:end-1)-x(3:end,3:end))/res2;
    
    neighfun{5,1} = @(x) (x(2:end-1,2:end-1)-x(3:end,2:end-1))/res;
    neighfun{6,1} = @(x) (x(2:end-1,2:end-1)-x(3:end,1:end-2))/res2;
    neighfun{7,1} = @(x) (x(2:end-1,2:end-1)-x(2:end-1,1:end-2))/res;
    neighfun{8,1} = @(x) (x(2:end-1,2:end-1)-x(1:end-2,1:end-2))/res2;
    
    % preallocate arrays
    G   = zeros(siz);
    
    
    
    % loop through neighbors
    if nargout==1;
        % if only gradient is required
        for neigh = (1:8);
            G        = max(neighfun{neigh,1}(DEM),G);
        end
    else
        ASP = G;
        % if gradient and aspect are required
        for neigh = (1:8);
            G2       = neighfun{neigh,1}(DEM);
            I        = G2>G;
            G(I)     = G2(I);
            ASP(I)   = neigh;
        end
    end
end

switch lower(unit)
    case 'tan'
        % do nothing
    case 'deg'
        G = atand(G);
    case 'rad'
        G = atan(G);
    case 'sin'
        G = sin(atan(G));
    otherwise
        error('unknown unit')
end
        


