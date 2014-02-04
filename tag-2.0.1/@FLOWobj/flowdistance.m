function OUT = flowdistance(FD,varargin) 

% flow distance in upstream and downstream direction
%
% Syntax
%
%     D = flowdistance(FD)
%     D = flowdistance(FD,GRIDobj)
%     D = flowdistance(FD,S)
%     D = flowdistance(FD,IX)
%     D = flowdistance(FD,x,y)
%     D = flowdistance(...,direction)
%
% Description
%
%     flowdistance calculates the horizontal distance along the flow
%     network in either upstream (default) or downstream direction.
%     Downstream flow distance is the maximum distance along the drainage
%     network (flow length). If only a flow direction object (FLOWobj) is
%     supplied to the function, flowdistance calculates the distance from
%     outlets and ridges in upstream and downstream direction,
%     respectively. If seed locations for the distance transform are
%     supplied (e.g. logical GRIDobj, STREAMobj, linear index or coordinate 
%     pairs) the distance is calculated from the seed locations. When using
%     a STREAMobj as input, usually upstream distance makes sense only.
%
% Input arguments
%
%     FD          flow direction (FLOWobj)
%     GRIDobj     logical grid (GRIDobj)
%     S           instance of STREAMobj
%     IX          linear index of seeds
%     x,y         x- and y-coordinate vectors
%     direction   'upstream' (default) or 'downstream' (maximum flow path
%                 length)
%
% Output arguments
%
%     D           distance grid (GRIDobj)
%
% Example
%
%
%
% See also: FLOWobj, FLOWobj/flowacc, GRIDobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 25. January, 2013


%% check input arguments

narginchk(1,4)

% check if last input argument is a string 
if ~isempty(varargin) && ischar(varargin{end})
    direction = validatestring(varargin{end},{'upstream','downstream'});
    nrargsin = nargin - 1;
else
    direction = 'upstream';
    nrargsin = nargin;
end



if nrargsin == 2;
    % SEED pixels are either supplied as logical matrix, GRIDobj, or linear
    % index
    SEED = varargin{1};
    isGRIDobj = isa(SEED,'GRIDobj');   
    if (islogical(SEED) || isGRIDobj)
        validatealignment(FD,SEED);
        if isGRIDobj;
            SEED = SEED.Z;
        end
    elseif isa(SEED,'STREAMobj');
        ix = SEED.IXgrid;
        SEED = false(FD.size);
        SEED(ix) = true;
        
        
    else
        % SEED is supposed to be supplied as linear index in the GRIDobj
        ix   = varargin{1};
        SEED = false(FD.size);
        SEED(ix) = true;
    end
elseif nrargsin == 3;
    % SEED pixels are supplied as coordinate pairs
    ix   = coord2ind(FD,varargin{1},varargin{2});
    SEED = false(FD.size);
    SEED(ix) = true;
end

if ~strcmpi(FD.type,'single')
    error('TopoToolbox:FLOWobj','flowdistance requires flow distance type to be single')
end

%% Do calculation
DIST = getdistance(FD.ix,FD.ixc,FD.size,FD.cellsize);

switch direction
    case 'upstream'
        %% Upstream distance calculation
        if nrargsin == 1;
            D     = zeros(FD.size);
            start = numel(FD.ix);
            for r = start:-1:1;
                D(FD.ix(r)) = D(FD.ixc(r))+DIST(r);
            end
        else
            D     = inf(FD.size);
            D(SEED) = 0;
            start = find(SEED(FD.ixc),1,'last');
            for r = start:-1:1;
                D(FD.ix(r)) = min(D(FD.ixc(r))+DIST(r),D(FD.ix(r)));
            end
            
            D(isinf(D)) = nan;
        end
    case 'downstream'
        
        %% Downstream distance calculation
        
        if nrargsin == 1
            D = zeros(FD.size);
        else
            D = -inf(FD.size);
            D(SEED) = 0;
        end
        
        for r = 1:numel(FD.ix);
            D(FD.ixc(r)) = max(D(FD.ix(r))+DIST(r),D(FD.ixc(r)));
        end
        
        if nrargsin >= 2;
            D(isinf(D)) = nan;
        end
end


%% Prepare Output
OUT = copy2GRIDobj(FD);
% write output to GRIDobj
OUT.Z = D;
OUT.zunit = '';
OUT.name  = [direction 'flow distance'];


