function DEMr = project(SOURCE,TARGET,varargin)

%PROJECT transforms a GRIDobj between projected coordinate systems
%
% Syntax
%
%     OUT = project(SOURCE,TARGET)
%     OUT = project(SOURCE,TARGET,pn,pv,...)
%
% Description
%
%     project reprojects the GRIDobj 'SOURCE' to have the same projection
%     as the GRIDobj 'TARGET'. Both GRIDobj's need to have an mstruct in
%     their 'georef' fields. If SOURCE does not, then it is assumed to have
%     a geographic coordinate system (horizontal WGS84 datum). The function
%     transforms the projected coordinates first to geographic coordinates
%     and then back to projected coordinates using the functions mfwdtran
%     and minvtran (both part of the Mapping Toolbox). The spatial
%     resolution of the output (OUT) will be the same as of TARGET, unless
%     set differently by the optional variable res.
%
% Input arguments
%
%     SOURCE         instance of GRIDobj that shall be transformed
%     TARGET         instance of GRIDobj with the target projection or 
%                    mapping projection structure (mstruct)
%
%     Parameter name/value pairs
%     
%     'res'          scalar
%                    spatial resolution. Default is TARGET.cellsize. If
%                    target is an mstruct, res must be supplied. 
%     'method'       string
%                    'bilinear' (default),'bicubic' or 'nearest'
%     'fillvalue'    scalar
%                    nan (for single and double grids). Otherwise 0.
%     'align'        logical scalar
%                    true (default) or false. If true, the OUT grid will be
%                    spatially aligned with TARGET. This means, they have
%                    the same extent and spatial alignment of cells. If
%                    'align', true then setting the resolution 'res' will
%                    have no effect.
%                    
% Output arguments
%
%     OUT            instance of GRIDobj
%
%
% See also: GRIDobj, reproject2utm, egm96heights
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de) and Wolfgang
%         Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 2. December, 2018


narginchk(2,inf)

if isa(TARGET,'GRIDobj')
    cs = TARGET.cellsize;
    targetisGRIDobj = true;
elseif isstruct(TARGET)
    % target is a mapping structure
    cs = [];
    targetisGRIDobj = false;
end

p = inputParser;
p.FunctionName = 'GRIDobj/project';
addParamValue(p,'res',  cs,@(x) isscalar(x));
addParamValue(p,'align',true,@(x) isscalar(x));
addParamValue(p,'method','bilinear',@(x) ischar(x));
addParamValue(p,'fillvalue',nan,@(x) isscalar(x));
parse(p,varargin{:});

% Get mstruct of SOURCE
try
    mstruct_s = SOURCE.georef.mstruct;
    if isempty(mstruct_s)
        error('error')
    end
    %cstype_s = SOURCEDEM.georef.SpatialRef.CoordinateSystemType;
    pr2pr     = true;
catch
    % there is either no mstruct available or the coordinate system is
    % unknown. We assume that the source coordinate system is geographic
    % (WGS84 Datum)
    pr2pr     = false;
end
    
% Get mstruct of TARGET
if targetisGRIDobj
    mstruct_t = TARGET.georef.mstruct;
    cs = p.Results.res;
else
    mstruct_t = TARGET;
    cs = p.Results.res;
    if isempty(cs)
        error('Spatial resolution of the output GRIDobj must be set.');
    end
end
    
    
    
%cstype_t = TARGETDEM.georef.SpatialRef.CoordinateSystemType;

if ~p.Results.align || ~targetisGRIDobj 

    [x0,y0]   = getcoordinates(SOURCE);
    x0        = x0(:);
    y0        = y0(:);

    % Calculate bounds of the reprojected DEM
    x         = [x0; x0; repmat(x0(1),numel(y0),1); repmat(x0(end),numel(y0),1)];
    y         = [repmat(y0(1),numel(x0),1); repmat(y0(end),numel(x0),1); y0; y0]; 
    
    if pr2pr
        [lat,lon] = minvtran(mstruct_s,x,y);
        [x,y]     = mfwdtran(mstruct_t,lat,lon);
    else
        [x,y]     = mfwdtran(mstruct_t,y,x);
    end
    
    tlims([1 2]) = [min(x),max(x)];
    tlims([3 4]) = [min(y),max(y)];
    
    res       = p.Results.res;
else
    
    [x0,y0]   = getcoordinates(SOURCE);
    [x,y]     = getcoordinates(TARGET);
    tlims([1 2]) = [min(x), max(x)];
    tlims([3 4]) = [min(y),max(y)];
    res       = TARGET.cellsize;
    
end

% Limits of source grid
slims([1 2])  = [min(x0) max(x0)];
slims([3 4])  = [min(y0) max(y0)];

% get fillvalue
fillval = p.Results.fillvalue;
if isinteger(SOURCE.Z)
    fillval = zeros(1,class(SOURCE.Z));
elseif islogical(SOURCE.Z)
    fillval = false;
else
    fillval = cast(fillval,class(SOURCE.Z));
end
fillval = double(fillval);
    
% check method
meth = validatestring(p.Results.method,{'bilinear','linear','nearest','bicubic'});

% Prepare tform for the image transform
if pr2pr
    T = maketform('custom', 2, 2, @FWDTRANSpr2pr, @INVTRANSpr2pr, []);
else
    T = maketform('custom', 2, 2, @FWDTRANSgcs2pr, @INVTRANSgcs2pr, []);
end

% Calculate image transform
[Znew,xdata,ydata] = imtransform(flipud(SOURCE.Z),T,meth,...
     'Xdata',tlims([1 2]),'Ydata',tlims([3 4]),...
     'Udata',slims([1 2]),'Vdata',slims([3 4]),...
     'XYScale',[res res],'Fillvalues',fillval);
 
Znew = flipud(Znew); 

if p.Results.align && targetisGRIDobj
    % If aligned, this is easy. The transformed GRIDobj will be perfectly
    % aligned with TARGET
    DEMr = TARGET;
    DEMr.Z = Znew;
else
    % We have calculated the imtransform with 'ColumnsStartFrom' south. 
    % GRIDobjs use 'ColumnsStartFrom' north

    xnew = cumsum([xdata(1) repmat(res,1,size(Znew,2)-1)]);
    ynew = flipud(cumsum([ydata(1) repmat(res,1,size(Znew,1)-1)])');

    % Construct GRIDobj
    DEMr = GRIDobj(xnew,ynew,Znew);
    % figure, subplot(1,2,1), imagesc(DEMr), subplot(1,2,2), imagesc(TARGETDEM)

    % Include geospatial and other information
    R = refmatToMapRasterReference(DEMr.refmat,DEMr.size);
    DEMr.georef.SpatialRef = R;
    DEMr.georef.mstruct = mstruct_t;
    if targetisGRIDobj
        DEMr.georef.GeoKeyDirectoryTag = TARGET.georef.GeoKeyDirectoryTag;
    else
        DEMr.georef.GeoKeyDirectoryTag = [];
        warning('Could not set GeoKeyDirectoryTag in output GRIDobj')
    end
    DEMr.zunit = SOURCE.zunit;
    DEMr.xyunit = SOURCE.xyunit;
end

DEMr.name = [SOURCE.name ' (projected)'];


% Transformation functions for imtransform (projected --> projected)
    function x = FWDTRANSpr2pr(u,~)
        [tlat,tlon] = minvtran(mstruct_s,u(:,1),u(:,2));
        [tx,ty] = mfwdtran(mstruct_t,tlat,tlon);
        x = [tx ty];        
    end

 
    function u = INVTRANSpr2pr(x,~)
        [tlat,tlon] = minvtran(mstruct_t,x(:,1),x(:,2));
        [tx,ty] = mfwdtran(mstruct_s,tlat,tlon);
        u = [tx ty];
    end

% Transformation functions for imtransform (gcs --> projected)
    function x = FWDTRANSgcs2pr(u,~)
        [tx,ty] = mfwdtran(mstruct_t,u(:,2),u(:,1));
        x = [tx ty];        
    end

 
    function u = INVTRANSgcs2pr(x,~)
        [tlat,tlon] = minvtran(mstruct_t,x(:,1),x(:,2));
        u = [tlon tlat];
    end

end