function c = idw(P,marks,varargin)

%IDW Inverse distance weighted interpolation on stream networks
%
% Syntax
%
%     c = idw(P,marks)
%     c = idw(P,marks,pn,pv,...)
%
% Description
%
%     idw computes an inverse distance weighted interpolation on a stream
%     network. Distances are calculated as geodesic distances on the
%     network (see graph/distances).
%
%     Note that idw is probably rarely a good choice as an interpolator on
%     a network. Better use STREAMobj/inpaintnans
%
% Input arguments
%
%     P       PPS object
%     marks   point attributes
%     
%     Parameter name/value pairs
%
%     'beta'       {2}. Exponent in the IDW equation
%     'chunksize'  Breaks computation into chunks so that the distance 
%                  matrix becomes not too large. A value of 1000 may be
%                  appropriate.
%     'extrap'     {false} or true. Don't extrapolate using idw. Hence, by
%                  default this option is set to false.
%
% Output arguments
%
%     c       node attribute list with interpolated values
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = removeshortstreams(S,100);
%     S = clean(S);
%     P = PPS(S,'rpois',0.001,'z',DEM);
%     marks = rand(npoints(P),1);
%     c = idw(P,marks);
%     plot(P.S,'k');
%     hold on
%     plotc(P,c)
%     hold off
% 
% See also: PPS, STREAMobj/inpaintnans, STREAMobj/crs, STREAMobj/smooth
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. March, 2020

distmatsize = 1000;

p = inputParser;
p.FunctionName = 'PPS/idw';
addRequired(p,'P',@(x) isa(x,'PPS'));
addRequired(p,'marks',@(x) numel(x) == npoints(P));

% Add elevation
addParameter(p,'beta',2,@(x) isscalar(x) && x > 0);
addParameter(p,'chunksize',round(distmatsize/npoints(P)*distmatsize));
addParameter(p,'extrap',false);

% Parse
parse(p,P,marks,varargin{:});

beta = p.Results.beta;
chunksize = p.Results.chunksize;
extrap = p.Results.extrap;

c = getnal(P.S);
c(:) = nan;
c(P.PP) = marks;

[CS,locb] = STREAMobj2cell(P.S);

% if there is only one drainage basin
if numel(CS) == 1
    c = networkidw(P.S,c,beta,chunksize);
else
    % if there are multiple drainage basins
    cc = cellfun(@(ix) c(ix),locb,'UniformOutput',false);
    for r = 1:numel(CS)
        cc{r} = networkidw(CS{r},cc{r},beta,chunksize);
    end
    
    % map back to actual nal
    for r = 1:numel(cc)
        c(locb{r}) = cc{r};
    end
end

if ~extrap
    I = isinf(netdist(P,'dir','up')) | isinf(netdist(P,'dir','down'));
    c(I) = nan;
end
    

end

function c = networkidw(S,c,beta,chunksize)
I = isnan(c);
tn = find(I);
sn = find(~I);

if isempty(sn)
    c(:) = nan;
    return
end

d = S.distance;
node_distance = abs(d(S.ix)-d(S.ixc));
G = graph(S.ix,S.ixc,node_distance);

n = numel(tn);
bins = ceil((1:n)'/chunksize);
ix = accumarray(bins,(1:numel(tn))',[],@(x) {tn(x)});
for r = 1:numel(ix)
    d = distances(G,sn,ix{r});
    d = d.^beta;
    w = 1./d;
    w = w./sum(w);
    c(ix{r}) = sum(c(sn).*w);
end
end

