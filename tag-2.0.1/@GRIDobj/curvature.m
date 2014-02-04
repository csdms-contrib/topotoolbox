function C = curvature(DEM,ctype)

% 8-connected neighborhood curvature of a digital elevation model 
%
% Syntax
%
%     C = curvature(DEM)
%     C = curvature(DEM,'type')
%
% Description
%     
%     curvature returns the second numerical derivative (curvature) of a
%     digital elevation model. By default, curvature returns the profile
%     curvature (profc). 
%
% Input arguments
%
%     DEM    digital elevation model (GRIDobj)
%     type   'profc' (default) : profile curvature 
%            'planc' : planform curvature
%            'tangc' : tangential curvature
%            'meanc' : mean curvature
%
% Output arguments
%
%     C      curvature (GRIDobj)
%
% Remarks
%     
%     Please note that curvature is not defined for cells with zero 
%     gradient. Here, curvature is set to zero.
%
%     All formulas are according to Olaya (2009) on page 151-152.
%
% Example
% 
%     
%
%
% Reference: Olaya, V. (2009): Basic land-surface parameters. In: 
%            Hengl, T. & Reuter, H. I. (Eds.), Geomorphometry. Concepts, 
%            Software, Applications, Elsevier, 33, 141-169.
%
% See also: GRADIENT8
%        
% Author:  Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. January, 2013


% check input arguments
narginchk(1,2);
if nargin == 1
    ctype = 'profc';
else
    ctype = validatestring(ctype,{'profc','planc','tangc','meanc'});
end

% create a copy of the DEM instance
C = DEM;
c = class(DEM.Z);
switch c
    case 'double'
        C.Z = double.empty(0,0);
    otherwise
        C.Z = single.empty(0,0);
end

% Large matrix support. Break calculations in chunks using blockproc
if numel(DEM.Z)>(5001*5001);
    blksiz = bestblk(size(DEM.Z),5000);    
    C.Z = blockproc(DEM.Z,blksiz,@(x) curvaturesub(x,C.cellsize),...
           'BorderSize',[1 1],...
           'Padmethod','symmetric',...
           'UseParallel',true);
else
    C.Z = curvaturesub(DEM.Z,C.cellsize);
end

C.name = ctype;


% subfunction

function curv = curvaturesub(dem,cs)

if isstruct(dem);
    dem = dem.data;
    % DEM has already been padded
    correctedges = false;
    shape = 'same';
else
    correctedges = true;
    shape = 'valid';
end
    
% First-order partial derivatives:
[p,q] = gradient(dem,cs);

if correctedges
    dem = padarray(dem,[1 1],'symmetric');
end
% Second order derivatives according to Evans method (see Olaya 2009)
%
% z1 z2 z3
% z4 z5 z6
% z7 z8 z9

% kernel for d2z/dx2
kernel = [1 -2 1; 1 -2 1; 1 -2 1]./(3*cs.^2);
r = conv2(dem,kernel,shape);
% kernel for d2z/dy2
kernel = kernel';
t = conv2(dem,kernel,shape);
% kernel for d2z/dxy
kernel = [-1 0 1; 0 0 0; 1 0 -1]./(4*cs.^2);
s = conv2(dem,kernel,shape);


%% Other options to calculate Second-order partial derivatives:
% r = gradient(p,cs);
% t = gradient(q',cs)';
% % Second-order mixed partial derivative:
% s = gradient(p',cs)';


switch ctype
    case 'profc'
        curv = - (p.^2 .* r + 2*p.*q.*s + q.^2.*t)./((p.^2 + q.^2).*(1 + p.^2 + q.^2).^(3));
    case 'tangc'
        curv = - (q.^2 .* r + 2*p.*q.*s + p.^2.*t)./((p.^2 + q.^2).*(1 + p.^2 + q.^2).^(3));
    case 'planc'
        curv = - (q.^2 .* r + 2*p.*q.*s + p.^2.*t)./((1 + p.^2 + q.^2).^(3));
    case 'meanc'
        curv = (p.^2 .* r + 2*p.*q.*s + q.^2.*t)./((p.^2 + q.^2).*(1 + p.^2 + q.^2)) ...
            - ((1+q).^2 .* r + 2*p.*q.*s + (1+p).^2.*t)./(2.*(1 + p.^2 + q.^2).^(3/2));
        
end

curv(isinf(curv) | isnan(curv)) = 0;
curv(isnan(dem)) = nan;
% curv = reshape(curv,size(dem));
% if correctedges
%     curv = curv(2:end-1,2:end-1);
% end

end
end