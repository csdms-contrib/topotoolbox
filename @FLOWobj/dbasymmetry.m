function [D,AS,ST] = dbasymmetry(FD,S,varargin)

%DBASYMMETRY Drainage basin asymmetry
%
% Syntax
%
%     [D,AS,ST] = dbasymmetry(FD,S)
%     [D,AS,ST] = dbasymmetry(FD,S,'pn',pv,...)
%
% Description
%
%     dbasymmetry calculates the drainage basin asymmetry. The function
%     derives the longest flow path within each drainage basin the outlet
%     of which is defined by the stream network S. It then calculates the 
%     part of each drainage basin that is hydrologically left to the  
%     longest flowpath. 
%
% Input arguments
%
%     FD      FLOWobj with single flow directions
%     S       STREAMobj 
%
%     Parameter name/values
%
%     'side'            {'right'} or 'left'
%     'extractlongest'  {true} or false. If false, dbasymmetry will take S
%                       to do the calculations. Thus, the output D will
%                       be the right or left contributing area to the
%                       stream network. If false, then the output ST will
%                       be the same as S.
%
% Output arguments
%
%     D       GRIDobj with ones on the hydrologically right side of the 
%             river network, 0.5 in river pixels, and zeros elsewhere (if
%             'side' is set to 'right')
%     AS      GRIDobj with asymmetry values for each drainage basin
%     ST      STREAMobj with longest flow path in each basin
%
% Example
%
%     DEM = readexample('taiwan');
%     DEM = inpaintnans(DEM);
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     [D,AS,ST] = dbasymmetry(FD,S);
%     % plot
%     subplot(1,2,1)
%     imageschs(DEM,D,'colorbar',false)
%     hold on
%     plot(ST,'w')
%     hold off
%     subplot(1,2,2)
%     imageschs(DEM,AS,'colormap',ttscm('vik'),'colorbar',true)
%     hold on
%     plot(ST,'w')
%     hold off
%  
% See also: STREAMobj, FLOWobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 6. December, 2021

p = inputParser;         
p.FunctionName = 'dbasymmetry';
addParameter(p,'side','right',@(x) ischar(validatestring(x,{'right','left'})));
addParameter(p,'extractlongest',true);
parse(p,varargin{:});

[DB,ixoutlet] = drainagebasins(FD,S);

if p.Results.extractlongest
    S   = extend2divide(S,FD);
    S   = trunk(S);

    % old version that failed to extract the longest stream reliably
    % DIST = flowdistance(FD,streampoi(S,'outlet','ix'));
    % IX  = 1:prod(FD.size);
    % I   = DB.Z > 0;
    
    % IX  = accumarray(DB.Z(I(:)),IX(I(:)),[],...
    %    @(ix) ix(find(DIST.Z(ix) == max(DIST.Z(ix)),1,'first')));
    % S   = STREAMobj(FD,'chan',IX);
end

if ~FD.fastindexing
    FD.fastindexing = true;
end

% Indexing madness to obtain two-edge segments
ix   = FD.ix;
ixc  = FD.ixc;
ixcc = ixc;
temp = FD.ixcix(ixc);
ixcc(temp ~= 0) = ixc(temp(temp ~= 0));

% coordinate matrices
[X,Y] = getcoordinates(GRIDobj(FD,'logical'),'matrix');

% calculate the angle for each double line segment with respect to ixc
% angle of 2 relative to 1= atan2(v2.y,v2.x) - atan2(v1.y,v1.x)
xr1 = (X(ix)  -X(ixc)); 
yr1 = (Y(ix)  -Y(ixc));
xr2 = (X(ixcc) - X(ixc));
yr2 = (Y(ixcc) - Y(ixc));

ang = atan2(yr1,xr1) - atan2(yr2,xr2);
ang(ang<0) = ang(ang <0) + 2*pi;

% Map angles to GRIDobj    
ANG = GRIDobj(FD,'single');
ANG.Z(ix) = ang;

% Extract angles at stream obj
s_ang = ANG.Z(S.IXgrid(S.ix));
ANG   = GRIDobj(FD,'single')*nan(1,'single');
% Write back to GRIDobj
ANG.Z(S.IXgrid(S.ixc)) = s_ang;

switch lower(p.Results.side)
    case 'right'
        I = ang <= ANG.Z(ixc);    
    case 'left'
        I = ang >= ANG.Z(ixc);
end

% cut flow network
ix(I)  = [];
ixc(I) = [];

% Get right or left drainage area
D = GRIDobj(FD,'single');
D.Z(S.IXgrid) = 1;
for r = numel(ixc):-1:1
    D.Z(ix(r)) = D.Z(ixc(r));
end
% Set stream pixels to .5
D.Z(S.IXgrid) = .5;

if nargout > 1
    % calculate flow accumulation
    A  = flowacc(FD);
    % calculate flow accumulation of right (or left) side only
    A2 = flowacc(FD,D);
    % calculate the ratio between both
    AS = A2/A;
    % extract values at outlets
    as = AS.Z(ixoutlet);
    % map values to entire drainage basins
    as = [nan;as];
    AS.Z = as(DB.Z+1);
    
    % return stream network
    ST = S;
end

