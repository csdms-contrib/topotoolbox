function [MS,varargout] = asymmetry(D,GRID)
% ASYMMETRY   directional asymmetry of divide segments
%
% Syntax
%
%     [MS] = asymmetry(D,DZ)
%     [MS,S] = asymmetry(D,DZ)
%
% Description
%
%     ASYMMETRY computes for each segment in the divide network the mean
%     hillslope relief asymmetry, represented as a vector that is centered
%     on the divide segment and points in the direction of lower hillslope
%     relief. Hillslope relief is computed from divide-adjacent pixels in 
%     the GRIDobj DZ, which should be calculated with the function 
%     'vertdistance2stream', and based on the same STREAMobj used for 
%     calculating D.
%     Note that the optional output 'S' contains values for all divide
%     nodes, as well as divide-segment averaged values distributed to all
%     divide nodes. It lends itself for diplaying using the function
%     'plotc'.
%
% Input
%
%     D         instance of class DIVIDEobj
%     DZ        instance of class GRIDobj, created with the function
%               'vertdistance2stream'
%
% Output
%
%     MS        mapping structure with POINT entries representing 
%               divide segments
%      .Geometry  - 'Point'
%      .X         - x coordinate 
%      .Y         - y coordinate 
%      .order     - order
%      .dist      - along-divide distance 
%      .u         - x-component of asymmetry (east is positive)
%      .v         - y-component of asymmetry (north is positive)
%      .theta     - angle from north of asymmetry direction
%      .rho       - magnitude of asymmetry
%
% Optional output     
%     
%     S         data structure with LINE entries representing divide
%               segments
%      .IX        - linear indices of divide segment nodes
%      .x         - x coordinate of divide nodes
%      .y         - y coordinate of divide nodes
%      .order     - order
%      .dist      - along-divide distance 
%      .u         - x-component of divide node asymmetry (east is positive)
%      .v         - y-component of divide node asymmetry (north is positive)
%      .theta     - angle from north of divide segment(!) asymmetry 
%      .rho       - magnitude of divide segment(!) asymmetry 
%
% Example
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,S);
%     D = divorder(D,'topo');
%     DZ = vertdistance2stream(FD,S,DEM);
%     [MS,S] = asymmetry(D,DZ);
%     for i = 1 : length(S)
%         S(i).length = max(getdistance(S(i).x,S(i).y));
%     end
%     imageschs(DEM,[],'colormap',[.9 .9 .9],'colorbar',false);
%     %imageschs(DEM,gradient8(DEM,'deg'),'caxis',[0 45]) % with slope as background
%     hold on
%     plotc(D,vertcat(S.rho),'caxis',[0 0.5],'limit',[1000 inf])
%     colormap(gca,flipud(pink))
%     axis image
%     hc = colorbar;
%     hc.Label.String = 'Divide asymmetry index';
%     hold on
%     ix = [MS.dist]>1000; % & [MS.rho]>0;
%     f = [S.length]./1e3;
%     quiver([MS(ix).X],[MS(ix).Y],[MS(ix).u].*f(ix),[MS(ix).v].*f(ix),2,...
%         'color','r','linewidth',1)
%     title('Drainage divide asymmetry and direction of lower hillslope relief')
%
%
% See also: DIVIDEobj, DIVIDEobj/sort
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: Nov 2018


% Get vectors 
[x,y] = ind2coord(D,vertcat(D.IX));

% Preprocess DZ grid
GRID.Z(GRID.Z<0) = 0;
GRID.Z(isinf(GRID.Z)) = nan;

% Calculate across divide attributes 
x1 = [NaN; x(1:end-1)];
x2 = [NaN; x(2:end)];
y1 = [NaN; y(1:end-1)];
y2 = [NaN; y(2:end)];
dx = x1-x2;
dy = y1-y2;
hcs = GRID.cellsize/2;
ix = dx==0; % vertical link
iy = dy==0; % horizontal link
meanx = (x1+x2)./2;
meany = (y1+y2)./2;
px = meanx + hcs.*ix;
qx = meanx - hcs.*ix; 
py = meany + hcs.*iy;
qy = meany - hcs.*iy;
pix = coord2ind(GRID,px,py); % top and right
qix = coord2ind(GRID,qx,qy); % bottom and left
nx = ~isnan(pix) & ~isnan(qix);

% direction of asymmetry (orthogonal to the orientation)
u = ix;
v = iy;

% magnitude of asymmetry
hr1 = nan(size(pix));
hr2 = hr1;
hr1(nx) = GRID.Z(pix(nx)); % top and right
hr2(nx) = GRID.Z(qix(nx)); % bottom and left
dhr = diff([hr1,hr2],1,2); % top minus bottom, right minus left
dhrn = dhr./sum([hr1,hr2],2); % 
u = u.*dhrn;
v = v.*dhrn;

S = onl2struct(D.IX,'x',x,'y',y,'order',D.order,'dist',D.distance,...
    'u',u,'v',v);
n = numel(S);

MS = struct('Geometry','Point',...
    'X',cell(n,1),...
    'Y',cell(n,1));

for i = 1 : length(MS)
    
    tx = S(i).x(1:end-1);
    ty = S(i).y(1:end-1);
    td = getdistance(tx,ty);
    if numel(td)>1
        MS(i).X = interp1(td,tx,max(td)/2);
        MS(i).Y = interp1(td,ty,max(td)/2);
    else
        MS(i).X = tx;
        MS(i).Y = ty;
    end
    MS(i).order = nanmean(S(i).order);
    MS(i).dist = nanmean(S(i).dist);
    MS(i).u = double(nanmean(S(i).u));
    MS(i).v = double(nanmean(S(i).v));
    
    [theta,rho] = cart2pol(MS(i).u,MS(i).v);
    theta = rad2deg(theta);
    theta = -theta+90;
    theta(theta<0) = theta(theta<0)+360;
    MS(i).theta = double(theta);
    MS(i).rho = double(rho);
    
    S(i).theta = [ones(size(tx)).*double(theta);nan];
    S(i).rho = [ones(size(tx)).*double(rho);nan];
    S(i).east = [ones(size(tx)).*double(MS(i).u);nan];
    S(i).north = [ones(size(tx)).*double(MS(i).v);nan];

end

if nargout>1
    varargout{1} = S;
end

