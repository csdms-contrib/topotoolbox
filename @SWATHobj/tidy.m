function OUT = tidy(SW,varargin)

% TIDYS removes overlapping points from SWATHobj
%
% Syntax
%
%     OUT = tidy(SW)
%     OUT = tidy(SW,FD)
%     OUT = tidy(SW,FD,'both')
%
% Description
%
%     TIDY(SW) uses the function 'lineSegmentIntersect' by U. Murat
%     Erdem, freely available on Matlab Central and included in this code,
%     to remove overlapping points in a SWATHobj by setting the
%     corresponding Z values to NaN.
%
%     TIDY(SW,FD) can be used when the SWATHobj SW was generated with a
%     STREAMobj derived from the FLOWobj FD and sets data points outside
%     the corresponding watershed to NaN.
%
%     TIDY(SW,FD,'both') does both of the above.
%
% Input arguments
%
%     SW     instance of SWATHobj
%     FD     instance of FLOWobj
%
% Output arguments
%
%     OUT    instance of SWATHobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     FA  = flowacc(FD);
%     S = STREAMobj(FD,FA>1000);
%     S = trunk(klargestconncomps(S,1));
%     SW = SWATHobj(DEM,S,'smooth',300,'plot',false);
%     SW = tidy(SW,FD,'both');
%     figure,plot(SW)
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: May, 2015


% check inputs
narginchk(1,3)

if ~isa(SW,'SWATHobj')
    error('First input needs to be of class SWATHobj')
end

if nargin>=2
    FD = varargin{1};
    if ~isa(FD,'FLOWobj')
        error('Second input needs to be of class FLOWobj')
    end
end
if nargin==3
    if ~strcmp('both',varargin{2})
        error('Unknown third input to function call')
    end
end


OUT = SW;

if nargin==1 || nargin==3
    
    % remove duplicates from mapping
    [ny,~] = size(SW.X);
    if mod(ny,2)==0; ctr = ny/2; xtr = 1;
    else ctr = (ny+1)/2; xtr = 0; end
    
    % do each side of the SWATHobj individually
    for sides = 1 : 2
        if sides==1; edge = 1;
        else edge = ny; ctr = ctr+1+xtr; end
        if edge < ctr; step = -1; else step = 1; end
        id = ctr:step:edge;
        d = abs(SW.disty(id));
        maxd = max(d);
        
        % the points near the center of the swath are defined as the starting points
        x1 = SW.X(ctr,:)';
        y1 = SW.Y(ctr,:)';
        % the points on the edge of the swath are defined as the ending points
        x2 = SW.X(edge,:)';
        y2 = SW.Y(edge,:)';
        XY1 = [x1,y1,x2,y2];
        XY1(isnan(XY1(:,1)),:) = [];
        XY2 = XY1;
        LSI = lineSegmentIntersect(XY1,XY2);
        U = triu(LSI.intAdjacencyMatrix);
        
        % check each intersection
        [nru,ncu] = size(U);
        IX = find(U);
        for k = 1 : length(IX)
            d1 = abs( LSI.intNormalizedDistance1To2(IX(k)) );
            d2 = abs( LSI.intNormalizedDistance2To1(IX(k)) );
            [p1,p2] = ind2sub([nru,ncu],IX(k));
            % determine which line to cut short
            if d1 <= d2;
                thisp = p2;
                thisd = d2*maxd;
            else
                thisp = p1;
                thisd = d1*maxd;
            end
            % exclude points along line beyond intersection
            ix = d>thisd;
            OUT.Z(id(ix),thisp) = nan;
            %OUT.X{k}(id(ix),thisp) = nan;
            %OUT.Y{k}(id(ix),thisp) = nan;
        end
    end

    
end



if exist('FD','var')
    
    FA = flowacc(FD);
    
    x = SW.xy0(:,1);
    y = SW.xy0(:,2);
    ix = coord2ind(FA,x,y);
    fa = FA.Z(ix);
    [~,ixfa] = sort(fa,'descend');
    ix = ix(ixfa);
    outl = ix(2);
    outl_fa = fa(2);
    
    [~,oix] = sort(outl_fa,'ascend');
    outl_s = outl(oix);
    
    ML = FLOWobj2GRIDobj(FD);
    ML.Z(:) = false;
    [x,y] = getcoordinates(ML);
    minx = min(x);
    maxx = max(x);
    miny = min(y);
    maxy = max(y);
    
    
    L = drainagebasins(FD,outl_s);
    
    ix_sw = 1:numel(OUT.X);
    x = OUT.X(ix_sw);
    y = OUT.Y(ix_sw);
    
    % out-of-range subscripts will be nan
    ix = x>=minx & x<=maxx & y>=miny & y<=maxy;
    ix_sw = ix_sw(ix);
    x = x(ix);
    y = y(ix);
    ix_grid = coord2ind(FD,x,y);
    ix_sw   = ix_sw(~isnan(ix_grid));
    ix_grid = ix_grid(~isnan(ix_grid));
    
    % set outside of watershed to nan
    nix = ~L.Z(ix_grid);
    OUT.Z(ix_sw(nix)) = nan;
    %OUT.Y{this_i}(ix_sw(nix)) = nan;
    %OUT.X{this_i}(ix_sw(nix)) = nan;
    
    nix = logical(ML.Z(ix_grid));
    OUT.Z(ix_sw(nix)) = nan;
    %OUT.Y{this_i}(ix_sw(nix)) = nan;
    %OUT.X{this_i}(ix_sw(nix)) = nan;
    


    
    
    
%     % collect outlets
%     outl = zeros(size(SW.Z));
%     outl_fa = zeros(size(SW.Z));
%     for i = 1:length(SW.Z)
%         x = SW.xy0{i}(:,1);
%         y = SW.xy0{i}(:,2);
%         ix = coord2ind(FA,x,y);
%         fa = FA.Z(ix);
%         [~,ixfa] = sort(fa,'descend');
%         ix = ix(ixfa);
%         outl(i) = ix(2);
%         outl_fa(i) = fa(2);
%     end
%     
%     [~,oix] = sort(outl_fa,'ascend');
%     outl_s = outl(oix);
% 
%     ML = FLOWobj2GRIDobj(FD);
%     ML.Z(:) = false;
%     for i = 1 : length(outl_s)
%         this_i = oix(i);
%         L = drainagebasins(FD,outl_s(i));
% 
%         ix_sw = 1:numel(OUT.X{this_i});
%         x = OUT.X{this_i}(ix_sw);
%         y = OUT.Y{this_i}(ix_sw);
% 
%         % set out-of-range subscripts to nan
%         xlims = FD.georef.SpatialRef.XLimWorld;
%         ylims = FD.georef.SpatialRef.YLimWorld;
%         ix = find(x>=xlims(1) & x<=xlims(2) & ...
%             y>=ylims(1) & y<=ylims(2));
%         ix_sw = ix_sw(ix);
%         x = x(ix);
%         y = y(ix);
% 
%         ix_grid = coord2ind(FD,x,y);
%         ix_sw   = ix_sw(~isnan(ix_grid));
%         ix_grid = ix_grid(~isnan(ix_grid));
% 
%         % set out-of-watershed to nan
%         nix = ~L.Z(ix_grid);
%         OUT.Z{this_i}(ix_sw(nix)) = nan;
%         %OUT.Y{this_i}(ix_sw(nix)) = nan;
%         %OUT.X{this_i}(ix_sw(nix)) = nan;
% 
%         nix = logical(ML.Z(ix_grid));
%         OUT.Z{this_i}(ix_sw(nix)) = nan;
%         %OUT.Y{this_i}(ix_sw(nix)) = nan;
%         %OUT.X{this_i}(ix_sw(nix)) = nan;
%         ML = ML|L;



    
    
    
end



end



%% Subfunctions
% The following code was downloaded from Matlab Central

function out = lineSegmentIntersect(XY1,XY2)
%LINESEGMENTINTERSECT Intersections of line segments.
%   OUT = LINESEGMENTINTERSECT(XY1,XY2) finds the 2D Cartesian Coordinates of
%   intersection points between the set of line segments given in XY1 and XY2.
%
%   XY1 and XY2 are N1x4 and N2x4 matrices. Rows correspond to line segments.
%   Each row is of the form [x1 y1 x2 y2] where (x1,y1) is the start point and
%   (x2,y2) is the end point of a line segment:
%
%                  Line Segment
%       o--------------------------------o
%       ^                                ^
%    (x1,y1)                          (x2,y2)
%
%   OUT is a structure with fields:
%
%   'intAdjacencyMatrix' : N1xN2 indicator matrix where the entry (i,j) is 1 if
%       line segments XY1(i,:) and XY2(j,:) intersect.
%
%   'intMatrixX' : N1xN2 matrix where the entry (i,j) is the X coordinate of the
%       intersection point between line segments XY1(i,:) and XY2(j,:).
%
%   'intMatrixY' : N1xN2 matrix where the entry (i,j) is the Y coordinate of the
%       intersection point between line segments XY1(i,:) and XY2(j,:).
%
%   'intNormalizedDistance1To2' : N1xN2 matrix where the (i,j) entry is the
%       normalized distance from the start point of line segment XY1(i,:) to the
%       intersection point with XY2(j,:).
%
%   'intNormalizedDistance2To1' : N1xN2 matrix where the (i,j) entry is the
%       normalized distance from the start point of line segment XY1(j,:) to the
%       intersection point with XY2(i,:).
%
%   'parAdjacencyMatrix' : N1xN2 indicator matrix where the (i,j) entry is 1 if
%       line segments XY1(i,:) and XY2(j,:) are parallel.
%
%   'coincAdjacencyMatrix' : N1xN2 indicator matrix where the (i,j) entry is 1
%       if line segments XY1(i,:) and XY2(j,:) are coincident.

% Version: 1.00, April 03, 2010
% Version: 1.10, April 10, 2010
% Author:  U. Murat Erdem

% CHANGELOG:
%
% Ver. 1.00:
%   -Initial release.
%
% Ver. 1.10:
%   - Changed the input parameters. Now the function accepts two sets of line
%   segments. The intersection analysis is done between these sets and not in
%   the same set.
%   - Changed and added fields of the output. Now the analysis provides more
%   information about the intersections and line segments.
%   - Performance tweaks.

% I opted not to call this 'curve intersect' because it would be misleading
% unless you accept that curves are pairwise linear constructs.
% I tried to put emphasis on speed by vectorizing the code as much as possible.
% There should still be enough room to optimize the code but I left those out
% for the sake of clarity.
% The math behind is given in:
%   http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
% If you really are interested in squeezing as much horse power as possible out
% of this code I would advise to remove the argument checks and tweak the
% creation of the OUT a little bit.

%%% Argument check.
%-------------------------------------------------------------------------------

validateattributes(XY1,{'numeric'},{'2d','finite'});
validateattributes(XY2,{'numeric'},{'2d','finite'});

[n_rows_1,n_cols_1] = size(XY1);
[n_rows_2,n_cols_2] = size(XY2);

if n_cols_1 ~= 4 || n_cols_2 ~= 4
    error('Arguments must be a Nx4 matrices.');
end

%%% Prepare matrices for vectorized computation of line intersection points.
%-------------------------------------------------------------------------------
X1 = repmat(XY1(:,1),1,n_rows_2);
X2 = repmat(XY1(:,3),1,n_rows_2);
Y1 = repmat(XY1(:,2),1,n_rows_2);
Y2 = repmat(XY1(:,4),1,n_rows_2);

XY2 = XY2';

X3 = repmat(XY2(1,:),n_rows_1,1);
X4 = repmat(XY2(3,:),n_rows_1,1);
Y3 = repmat(XY2(2,:),n_rows_1,1);
Y4 = repmat(XY2(4,:),n_rows_1,1);

X4_X3 = (X4-X3);
Y1_Y3 = (Y1-Y3);
Y4_Y3 = (Y4-Y3);
X1_X3 = (X1-X3);
X2_X1 = (X2-X1);
Y2_Y1 = (Y2-Y1);

numerator_a = X4_X3 .* Y1_Y3 - Y4_Y3 .* X1_X3;
numerator_b = X2_X1 .* Y1_Y3 - Y2_Y1 .* X1_X3;
denominator = Y4_Y3 .* X2_X1 - X4_X3 .* Y2_Y1;

u_a = numerator_a ./ denominator;
u_b = numerator_b ./ denominator;

% Find the adjacency matrix A of intersecting lines.
INT_X = X1+X2_X1.*u_a;
INT_Y = Y1+Y2_Y1.*u_a;
INT_B = (u_a >= 0) & (u_a <= 1) & (u_b >= 0) & (u_b <= 1);
PAR_B = denominator == 0;
COINC_B = (numerator_a == 0 & numerator_b == 0 & PAR_B);


% Arrange output.
out.intAdjacencyMatrix = INT_B;
out.intMatrixX = INT_X .* INT_B;
out.intMatrixY = INT_Y .* INT_B;
out.intNormalizedDistance1To2 = u_a;
out.intNormalizedDistance2To1 = u_b;
out.parAdjacencyMatrix = PAR_B;
out.coincAdjacencyMatrix= COINC_B;

end

