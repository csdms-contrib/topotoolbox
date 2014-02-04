function OUT = dependencemap(FD,varargin)


% upslope area for specific locations in a digital elevation model
%
% Syntax
%
%     I = dependencemap(FD,L)
%     I = dependencemap(FD,ix)
%     I = dependencemap(FD,x,y)
%
% Description
%
%     dependencemap returns a logical matrix masking the upslope part
%     of the digital elevation model that drain into the specified cells.
%
% Input
%
%     FD        flow direction object (FlowDirObj)
%     L         logical grid (GRIDobj)
%     ix        linear index into a grid
%     x,y       coordinate pairs
%
% Output
%
%     I         logical influence grid (GRIDobj)
%
% Example
%
%
% See also: FLOWobj, FLOWobj/INFLUENCEMAP
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. January, 2013



%% check input arguments
narginchk(1,3)
if nargin == 2;
    % SEED pixels are either supplied as logical matrix, GRIDobj, or linear
    % index
    SEED = varargin{1};
    isGRIDobj = isa(SEED,'GRIDobj');   
    if (islogical(SEED) || isGRIDobj)
        validatealignment(FD,SEED);
        if isGRIDobj;
            SEED = SEED.Z;
        end
    else
        % SEED is supposed to be supplied as linear index in the GRIDobj
        ix   = varargin{1};
        SEED = false(FD.size);
        ix   = round(ix);
        if any(ix <= 0 | ix > prod(FD.size));
            error('TopoToolbox:WrongInput',...
            ['Linear indices must not be less or equal to zero \n' ...
             'or larger than ' double2str(prod(FD.size)) '.']);
        end
        
        SEED(ix) = true;
    end
elseif nargin == 3;
    % SEED pixels are supplied as coordinate pairs
    ix   = coord2ind(FD,varargin{1},varargin{2});
    SEED = false(FD.size);
    SEED(ix) = true;
end


% this is a very crude and slow implementation of a graph traversal
% algorithm since all nodes are visited

if ~(exist(['dependencemap_mex.' mexext],'file') == 3);
    % m implementation
    for r = numel(FD.ix):-1:1;
        SEED(FD.ix(r)) = SEED(FD.ix(r)) || SEED(FD.ixc(r));
    end
    
else
    % mex implementation
    SEED = dependencemap_mex(FD.ix,FD.ixc,FD.size,find(SEED));    
end


%% Prepare Output
% empty GRIDobj
OUT = copy2GRIDobj(FD);
% write output to GRIDobj
OUT.Z = SEED;
OUT.zunit = 'logical';
OUT.name  = 'dependence map';


end








