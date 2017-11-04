function [DEMr] = projectGRIDobj(SOURCE,TARGET,res)

%PROJECTGRIDOBJ reprojects a GRIDobj
%
% Syntax
%
%     OUT = projectGRIDobj(SOURCE,TARGET)
%     OUT = projectGRIDobj(SOURCE,TARGET,res)
%
% Description
%
%     projectGRIDobj reprojects the GRIDobj 'SOURCE' to have the same
%     projection as the GRIDobj 'TARGET'. Both GRIDobj's need to have an
%     mstruct in their 'georef' fields. The function transforms the
%     projected coordinates first to geographic coordinates and then back
%     to projected coordinates using the functions mfwdtran and minvtran
%     (both part of the Mapping Toolbox). The spatial resolution of the 
%     output (OUT) will be the same as of SOURCE, unless set differently by
%     the optional variable res.
%
% Input arguments
%
%     SOURCE         instance of GRIDobj that shall be transformed
%     TARGET         instance of GRIDobj with the target projection
%     res            image resolution in map units
% 
% Output arguments
%
%     OUT            instance of GRIDobj
%
%
% See also: GRIDobj, reproject2utm
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 22. August, 2017


narginchk(2,3)

if nargin==2
    res = SOURCE.cellsize;
end


% Get mstruct of SOURCDEM
mstruct_s = SOURCE.georef.mstruct;
%cstype_s = SOURCEDEM.georef.SpatialRef.CoordinateSystemType;
 
% Get mstruct of TARGETDEM
mstruct_t = TARGET.georef.mstruct;
%cstype_t = TARGETDEM.georef.SpatialRef.CoordinateSystemType;


% Calculate bounds of the reprojected DEM
[x0,y0] = getcoordinates(SOURCE);
[X,Y] = meshgrid(x0,y0);
[lat,lon] = minvtran(mstruct_s,X(:),Y(:));
[x,y] = mfwdtran(mstruct_t,lat,lon);
lims([1 2]) = [min(x),max(x)];
lims([3 4]) = [min(y),max(y)];
 
% Prepare tform for the image transform
T = maketform('custom', 2, 2, @FWDTRANS, @INVTRANS, []);

% Calculate image transform
[Znew,xdata,ydata] = imtransform(flipud(SOURCE.Z),T,'linear',...
     'Xdata',lims([1 2]),'Ydata',lims([3 4]),...
     'Udata',x0([1 end]),'Vdata',y0([end 1])',...
     'XYScale',[res res],'Fillvalues',nan);
 
% We have calculated the imtransform with 'ColumnsStartFrom' south. 
% GRIDobjs use 'ColumnsStartFrom' north
Znew = flipud(Znew);
xnew = cumsum([xdata(1) repmat(res,1,size(Znew,2)-1)]);
ynew = flipud(cumsum([ydata(1) repmat(res,1,size(Znew,1)-1)])');

% Construct GRIDobj
DEMr = GRIDobj(xnew,ynew,Znew);
% figure, subplot(1,2,1), imagesc(DEMr), subplot(1,2,2), imagesc(TARGETDEM)

% Include geospatial and other information
R = refmatToMapRasterReference(DEMr.refmat,DEMr.size);
DEMr.georef.SpatialRef = R;
DEMr.georef.mstruct = TARGET.georef.mstruct;
DEMr.georef.GeoKeyDirectoryTag = TARGET.georef.GeoKeyDirectoryTag;
DEMr.zunit = SOURCE.zunit;
DEMr.xyunit = SOURCE.xyunit;
DEMr.name = [SOURCE.name ' (projected)'];


% Transformation functions for imtransform
    function x = FWDTRANS(u,~)
        [tlat,tlon] = minvtran(mstruct_s,u(:,1),u(:,2));
        [tx,ty] = mfwdtran(mstruct_t,tlat,tlon);
        x = [tx ty];        
    end

 
    function u = INVTRANS(x,~)
        [tlat,tlon] = minvtran(mstruct_t,x(:,1),x(:,2));
        [tx,ty] = mfwdtran(mstruct_s,tlat,tlon);
        u = [tx ty];
    end
end