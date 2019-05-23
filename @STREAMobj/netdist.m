function d = netdist(S,a,varargin)

%NETDIST distance transform on a stream network
%
% Syntax
%
%     d = netdist(S,a)
%     d = netdist(S,ix)
%     d = netdist(S,I)
%     d = netdist(S,c)
%     d = netdist(...,pn,pv,...)
% 
% Description
%
%     netdist computes the distance along a stream network S. It calculates
%     the distance of each node in the network from the nearest non-zero
%     node in the node-attribute list a.
%
% Input arguments
%
%     S     STREAMobj
%     a     logical node-attribute list
%     ix    index into GRIDobj
%     I     logical GRIDobj
%     c     all strings accepted by the function streampoi. For example,
%           d  = netdist(S,'outlet','dir','up') is the same as S.distance.
%
%     Parameter name/value pairs
%
%     'dir'     'both' (default) calculates the distance in upstream and
%               downstream direction. 'up' calculates distances only in
%               upstream direction and 'down' in downstream direction.
% 
%     'split'   false (default) or true. If true, netdist will do the
%               computation in parallel on each drainage basin. If the
%               Parallel Computing Toolbox is available, that may be 
%               faster. However, splitting the network in its connected
%               components produces some considerable computational
%               overhead.
%
% Output arguments
%
%     d     node-attribute list with distances. Nodes that cannot be
%           reached will have a distance of inf. 
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',500);
%     IX  = randlocs(S,100);
%     [x,y] = ind2coord(DEM,IX);
%     d   = netdist(S,IX);
%     plotc(S,d)
%     hold on
%     plot(x,y,'ok');
%
%
%
% See also: STREAMobj/distance, 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. May, 2019

% Check input arguments
p = inputParser;
p.FunctionName = 'STREAMobj/netdist';
addParameter(p,'split',false,@(x) isscalar(x));
addParameter(p,'dir','both');
parse(p,varargin{:});

% validate direction
direction = validatestring(p.Results.dir,{'both','up','down'});

% handle input
if isa(a,'GRIDobj')
    a = logical(getnal(S,a));
elseif ischar(a)
    a = streampoi(S,a,'logical');
elseif ~isnal(S,a)
    a = ismember(S.IXgrid,a);
elseif isnal(S,a)
    a = logical(a);
end

% run in parallel
if ~p.Results.split
    % computations are done in netdistsub
    d = netdistsub(S,a,direction);
else
    [CS,locb] = STREAMobj2cell(S);
    nrbasins = numel(CS);
    if nrbasins == 1
        d = netdistsub(S,a,direction);
    else
        D = cell(nrbasins,1);
        A = cellfun(@(ix) a(ix),locb,'UniformOutput',false);
        parfor r = 1:nrbasins
            D{r} = netdistsub(CS{r},A{r},direction);
        end
        d = zeros(size(a));
        for r = 1:nrbasins
            d(locb{r}) = D{r};
        end
        
    end
end
end

function d = netdistsub(S,a,direction)

ix  = S.ix;
ixc = S.ixc;
dd  = hypot(S.x(ix)-S.x(ixc),S.y(ix)-S.y(ixc));
d   = inf(size(a));
d(logical(a)) = 0;

switch direction
    case {'down','both'}
        for r = 1:numel(S.ix)
            d(ixc(r)) = min(d(ix(r))+dd(r),d(ixc(r)));
        end
end

switch direction
    case {'up','both'}
        
        for r = numel(ix):-1:1
            d(ix(r)) = min(d(ixc(r))+dd(r),d(ix(r)));
        end
end
    
end