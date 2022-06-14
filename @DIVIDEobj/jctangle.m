function [J,varargout] = jctangle(D,varargin)
%JCTANGLE   angles between divide segments at junctions
%
% Syntax
%
%     J = jctangle(D)
%     J = jctangle(D,unit)
%     [J,L] = jctangle(D,...)
%
% Description
%
%     Junctions in a divide network are points, where two or 
%     more divide segments meet. JCTANGLE measures the angle
%     between adjacent divide segments. In the current version,
%     only junctions of three divide segments are considered,
%     which are the most common, anyway. 
%
% Input
%
%     D        instance of class DIVIDEobj
%     unit     angular unit {'deg','rad'}
%
% Output
%
%     J    structure with the following fields:
%      .IX     linear index of junction
%      .tag    flag if angle data for this junction exists
%      .x,y    coordinates of junction
%      .a,b,c  angles between divide segments
%      .t1ix,t2ix,t3ix  linear index of joining divide segments
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,ST);
%     D = divorder(D,'topo');
%     J = jctangle(D);
%     A = [[J.a]',[J.b]',[J.c]'];
%     maxa = max(A,[],2);
%     plot(D,'color',[.8 .8 .8],'limit',[0 inf])
%     hold on
%     scatter([J.x],[J.y],80,maxa,'filled')
%     hold off
%     axis image
%     hc = colorbar;
%     hc.Label.String = 'Maximum junction angle (deg)';
%
% See also: DIVIDEobj, DIVIDEobj/sort
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020


if nargin==2
    unit = lower(varargin{1});
else
    unit = 'd';
end

switch unit
    case 'degree';  unit = 'deg';
    case 'deg';     unit = 'deg';
    case 'd';       unit = 'deg';
    case 'radians'; unit = 'rad';
    case 'rad';     unit = 'rad';
    case 'r';       unit = 'rad';
end


ix = D.jctedg==3;
ixjct = D.jct(ix);
[x,y] = ind2coord(D,ixjct);

ix = num2cell(ixjct);
x = num2cell(x);
y = num2cell(y);


J = struct('IX',ix,'tag',false,'x',x,'y',y,...
    'a',nan,'b',nan,'c',nan,...
    't1ix',nan,'t2ix',nan,'t3ix',nan);
%     't1d',nan,'t2d',nan,'t3d',nan);

M = onl2struct(D.IX);

% Loop over junctions
for i = 1 : length(ixjct)
    
    tjct = ixjct(i);
    % Find divide segments connected to the junction
    [ir,ic] = find(vertcat(M.st)==tjct);
    
    % Get divide segments and flip if necessary
    ix1 = M(ir(1)).IX(1:end-1);
    if ic(1)==2; ix1=flipud(ix1); end
    ix2 = M(ir(2)).IX(1:end-1);
    if ic(2)==2; ix2=flipud(ix2); end
    ix3 = M(ir(3)).IX(1:end-1);
    if ic(3)==2; ix3=flipud(ix3); end
    
    % Convert to coordinates and get length
    [tx,ty] = ind2coord(D,ix1);
	dist1 = cumsum([0;hypot(diff(tx),diff(ty))]);
    [tx,ty] = ind2coord(D,ix2);
    dist2 = cumsum([0;hypot(diff(tx),diff(ty))]);
    [tx,ty] = ind2coord(D,ix3);
    dist3 = cumsum([0;hypot(diff(tx),diff(ty))]);

    J(i).t1ix = ix1;
    J(i).t2ix = ix2;
    J(i).t3ix = ix3;
%     J(i).t1d = dist1;
%     J(i).t2d = dist2;
%     J(i).t3d = dist3;
    J(i).tag = true;
end


% Get angle data
[x,y] = refmat2XY(D.refmat,D.size);
[X,Y] = meshgrid(x,y);

L = struct; ct = 0;
for i = 1 : length(J)
    if J(i).tag
        
        ix = J(i).t1ix;
        tx = mean(X(ix)-X(ix(1)));
        ty = mean(Y(ix)-Y(ix(1)));
        [theta1,~] = cart2pol(tx,ty);
        rho1 = hypot(diff(X(ix([1,end]))),diff(Y(ix([1,end]))));
        
        ix = J(i).t2ix;
        tx = mean(X(ix)-X(ix(1)));
        ty = mean(Y(ix)-Y(ix(1)));
        [theta2,~] = cart2pol(tx,ty);
        rho2 = hypot(diff(X(ix([1,end]))),diff(Y(ix([1,end]))));
        
        ix = J(i).t3ix;
        tx = mean(X(ix)-X(ix(1)));
        ty = mean(Y(ix)-Y(ix(1)));
        [theta3,~] = cart2pol(tx,ty);
        rho3 = hypot(diff(X(ix([1,end]))),diff(Y(ix([1,end]))));
        
        % Sort and take difference
        t = [theta1,theta2,theta3];
        st = sort(t);
        dst = diff(st);
        J(i).a = dst(1);
        J(i).b = dst(2);
        J(i).c = 2*pi-sum(dst);
        
        % Lines that indicate the average orientation
        ct = ct+1;
        [dx,dy] = pol2cart(theta1,rho1);
        L(ct).x1 = [J(i).x;J(i).x+dx;nan];
        L(ct).y1 = [J(i).y;J(i).y+dy;nan];
        [dx,dy] = pol2cart(theta2,rho2);
        L(ct).x2 = [J(i).x;J(i).x+dx;nan];
        L(ct).y2 = [J(i).y;J(i).y+dy;nan];
        [dx,dy] = pol2cart(theta3,rho3);
        L(ct).x3 = [J(i).x;J(i).x+dx;nan];
        L(ct).y3 = [J(i).y;J(i).y+dy;nan];
        
    end
end



if strcmp(unit,'deg')
    a = num2cell(rad2deg([J.a]));
    [J.a] = a{:};
    b = num2cell(rad2deg([J.b]));
    [J.b] = b{:};
    c = num2cell(rad2deg([J.c]));
    [J.c] = c{:};
end

if nargout>1 && isfield(L,'x1')
    clear M
    M.x1 = vertcat(L.x1);
    M.y1 = vertcat(L.y1);
    M.x2 = vertcat(L.x2);
    M.y2 = vertcat(L.y2);
    M.x3 = vertcat(L.x3);
    M.y3 = vertcat(L.y3);
    varargout{1} = M;
else
    varargout{1} = [];
end


end

