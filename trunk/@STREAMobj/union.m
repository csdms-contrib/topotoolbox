function S = union(varargin)

% merge different instances of STREAMobj into a new instance
%
% Syntax
%
%     S = union(S1,S2,...,FD)
%
% Description
%
%     union combines different instances of STREAMobj into a new STREAMobj.
%     Note that this requires an instance of FLOWobj as last input
%     variable and that it is assumed that all STREAMobjs have been derived
%     from this FLOWobj (e.g., all STREAMobjs are subgraphs of the
%     FLOWobj). If this is not the case, the function might have an
%     unexpected behavior. In terms of set theory, the functions produces
%     a union of the nodes in all instances of S1.
%
% Input arguments
%
%     S1, S2, ...    several instances of STREAMobj
%     FD             instance of FLOWobj
%
% Output arguments
%
%     S              STREAMobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     A  = flowacc(FD);
%     S  = STREAMobj(FD,A>1000);
%     CS = STREAMobj2cell(S);
%     Sn = union(CS{1:2:end},FD);
%     plot(S,'r');
%     hold on
%     plot(Sn,'b');
%
% See also: STREAMobj, FLOWobj, STREAMobj/modify, STREAMobj/trunk
%           STREAMobj/STREAMobj2cell, STREAMobj/intersect
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. October, 2013

narginchk(2,inf)

% are instances spatially aligned?
for r = 1:numel(varargin);
    validatealignment(varargin{1},varargin{r})
end

% is the last input argument a FLOWobj
if ~isa(varargin{end},'FLOWobj')
    error('TopoToolbox:wronginput',...
        'The last input argument must be an instance of FLOWobj');
end

% empty GRIDobj
G = GRIDobj([]);
% find common properties of F and G and from F to G
pg = properties(G);
pf = properties(varargin{end});
p  = intersect(pg,pf);
for r = 1:numel(p);
    G.(p{r}) = varargin{end}.(p{r});
end
G.Z = false(G.size);

% create a stream raster
for r = 1:numel(varargin)-1;
    G.Z(varargin{r}.IXgrid) = true;
end

% use stream raster to create a new STREAMobj
S = STREAMobj(varargin{end},G);


    


