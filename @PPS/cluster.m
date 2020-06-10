function c = cluster(P,varargin)

%CLUSTER Hierarchical clustering of points in PPS
%
% Syntax
%
%     c = cluster(P)
%     c = cluster(P,pn,pv,...)
%
% Description
%
%     cluster uses hierarchical clustering to spatially classify the points
%     in the PPS-object P based on the interpoint-distance along the stream
%     network in P.
%
%     Note that cluster forms the full distance matrix which puts
%     constraints on the maximum amount of points in P. 
%
% Input arguments
%
%     P    instance of PPS
%
%     Parameter name value pairs
%
%     'type'        'interpoint' (default) or 'marks'
%     'cutoff'      threshold for creating clusters (the unit is according
%                   to the distance metric applied)
%
%     if 'type' is 'interpoint', following parameter name values apply
%
%     'val'         Node attribute list of distance values, if distance
%                   should be based on a different metric. Default is
%                   P.S.distance.
%     'd3d'         If true, distances are calculated in 3d. Default is
%                   false.
%
%     if 'type' is 'marks'
%
%     'val'         n*m matrix with values for each point (n=npoints(P)). 
%     'distance'    distance metric for calculating dissimilarities between
%                   mark values (see pdist)
%
%     parameter name/value pairs for both are
%
%     'method'      see function linkage. Default is 'average'
%     'criterion'   see function cluster. Default is 'distance'
%     'maxclust'    if set, than the function ignores the parameters
%                   'criterion' and 'cutoff'. See function cluster.
%
% Output arguments
%
%     c    cluster assignments for each point in P
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     % create six random points which serve as spatial cluster centers
%     rng(1); % make results replicable
%     P = PPS(S,'runif',6,'z',DEM);
%     d = density(P,'bw',400);
%     P = simulate(P,'intensity',d*10);
%     c = cluster(P,'maxclust',6);
%     convhull(P,'group',c)
%
%
% See also: PPS, PPS/convhull
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 11. February, 2019

% Check input arguments
p = inputParser;
p.FunctionName = 'PPS/cluster';
addRequired(p,'P',@(x) npoints(x) > 0)
addParameter(p,'type','interpoint')
addParameter(p,'method','average');
addParameter(p,'distance','euclidean');
addParameter(p,'DistParameter',[])
addParameter(p,'criterion','distance');
addParameter(p,'cutoff',tlength(P)/npoints(P));
addParameter(p,'maxclust',[],@(x) isscalar(x) & x > 0 );
addParameter(p,'val',[]);
addParameter(p,'d3d',false);
% Parse
parse(p,P,varargin{:});

type = validatestring(p.Results.type,{'interpoint','marks'});

% Calculate distances

switch type
    case 'interpoint'
        d = pointdistances(P,'output','matrix',...
            'type','graph',...
            'd3d',p.Results.d3d,...
            'val',p.Results.val);
        
        d = squareform(d,'tovector');
    case 'marks'
        if size(p.Results.val,1) ~= npoints(P)
            error('PPS:cluster','Wrong size of input parameter ''val''')
        end
        d = pdist(p.Results.val,p.Results.distance,p.Results.DistParameter);
        
end

% Agglomerative hierarchical cluster tree
L = linkage(d,p.Results.method);

% Extract clusters
if isempty(p.Results.maxclust)
    c = cluster(L,'cutoff',p.Results.cutoff,'criterion','distance');
else
    c = cluster(L,'maxclust',p.Results.maxclust);
end

end

