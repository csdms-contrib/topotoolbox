function M = flowdir(DEM,varargin)

% multiple and single flow direction algorithm for Digital Elevation Models
%
% Syntax
%
%     M = flowdir(DEM)
%     M = flowdir(DEM,'propertyname',propertyvalue,...)
%     [M,W] = flowdir(...)
%
% Description
% 
%     Multiple and single flowdirection algorithm that routes 
%     through flat terrain (not sinks).
%
% Input
%
%     X,Y       coordinate matrices created by meshgrid
%     dem       digital elevation model same size as X and Y
%
% Properties
%
% propertyname     propertyvalues
%
%     'type'            'multi' (default): multiple flowdirection (dinf)
%                       'single': single flow direction (d8). Flow occurs 
%                       only along the steepest descent
% 
%     'exponent'        exponent governing the relation between flow
%                       direction and slope. Default is 1.1, which means,  
%                       there is a nonlinear relation. You may want to 
%                       increase the exponent when flow direction should 
%                       rather follow a steepest descent (single) flow 
%                       direction (e.g. 5). This option is only effective 
%                       for multiple flowdirection.
%
%     'routeflats'      choose method to route over flats/plateaus.
%                       'route' uses the function routeflats
%                       'geodesic' uses routegeodesic (requires Matlab
%                       2011b or higher)
%                       'none' does not apply any routing through flats
%
% Output
%
%     M         flowdirection (sparse matrix)
%
% See also: FLOWobj
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 15. March, 2009



p = inputParser;
p.FunctionName = 'flowdir';


addParameter(p,'exponent',1.1,@(x) isscalar(x));
addParameter(p,'routeflats','geodesic'); %{'route','geodesic','none'};
addParameter(p,'type','multi');

parse(p,varargin{:});


% *********************************************************************
% normal multiple FLOW DIRECTION calculation

% calculate maximum slope and slope direction
% find neighbors of cells
[ic,icd] = ixneighbors(DEM.Z);
I = (DEM.Z(icd)-DEM.Z(ic)) > 0;
ic(I) = [];
icd(I) = [];
e = (DEM.Z(ic)-DEM.Z(icd))./getdistance(ic,icd,DEM.size,DEM.cellsize);
nrc = numel(DEM.Z);
% *********************************************************************
% flow direction matrix
M = sparse(ic,icd,double(e),nrc,nrc);


% *********************************************************************
% routing through flats
%
% A flat or plateau exists when two or more adjacent cells have the same 
% elevation. Up to now flow direction indicates for these cells
% no exchange of water.
% The subsequent code first identifies flats and then iteratively removes
% them.

switch p.Results.routeflats
    case 'route'
        [icf,icn] = routeflats(DEM.Z,p.Results.type);
        M = sparse(icf,icn,1,nrc,nrc)+M;        
    case 'geodesic'
        [icf,icn] = routegeodesic(DEM,p.Results.type);
        M = sparse(icf,icn,1,nrc,nrc)+M;
    case 'none'
end

% ******************************************************************
% single flow direction, flow concentration
switch p.Results.type
    case 'single'
        [m,IX2] = max(M,[],2);
        i = m==0;
        IX1 = (1:nrc)';
        IX1(i) = [];
        IX2(i) = [];
        M = sparse(IX1,IX2,1,nrc,nrc);
    otherwise
        if p.Results.exponent ~= 1;
            M = spfun(@(x) x.^p.Results.exponent,M);
        end
end

% ******************************************************************
% Row standardization of M only necessary when multiple flow dir
switch p.Results.type
    case 'multi'
        M = spdiags(spfun(@(x) 1./x,sum(M,2)),0,nrc,nrc) * M;
end

% and this is it...



