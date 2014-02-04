function [varargout] = baselevel(SW,varargin)
% create a SWATHobj with the local across-swath minimum
%
% Syntax
%
%     OUT = baselevel(SW)
%     OUT = baselevel(SW,DEM)
%     OUT = baselevel(SW,DEM,'pn','pv',...)
%     [OUT,DZ] = baselevel(SW,...)
%
% Description
%
%     baselevel(SW) creates a SWATHobj the same extent as SW and populates
%     the z-value fields with the across-swath minimum as determined from
%     the z-values in the field 'SW.zd0'.
%
%     baselevel(SW,DEM) creates a GRIDobj the same extent as DEM, where
%     pixels covered by the SWATHobj SW are given the across-swath minimum
%     z-values as obtained from DEM at the positions given in 'SW.xy0'.
%
% Input arguments
%
%     SW     instance of SWATHobj
%     DEM    instance of GRIDobj
%
%     Parameter name/value pairs
%
%     'values'    {[]}, n x 1 or 2 array
%     this parameter can be used to provide z-values which are used as the
%     local minimum. The values have to come as the output of the function
%     STREAMobj2XY. That means, the items of the SWATHobj separated by 
%     NaN's and ordered the same way as in the SWATHobj, and the values 
%     from each item ordered with increasing distance along the SWATHobj.
%     Distance values can also be provided and have to follow in the second
%     column of the array.
%
%     'smooth'    {0}, scalar
%     filter length in map units, which is used to smooth the z-values
%     along the direction of the SWATHobj using Matlab's filtfilt function
%
% Output arguments
%
%     OUT    instance of SWATHobj or GRIDobj
%     DZ     n x 2 array with distance-elevatio pairs for the line segments
%            of the SWATHobj. Individual segments are separated by 'NaNs'
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,A>1000);
%     S = klargestconncomps(S,1);
%     S = removeshortstreams(S,1e3);
%     SW = SWATHobj(DEM,S,'smooth',300,'plot',false);
%     SW = tidyswath(SW,FD,'both');
%     BL = baselevel(SW,DEM);
%     figure,imagesc(DEM-BL)
%     title('Elevation above river level (m)')
%
% See also: SWATHobj2GRIDobj
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: November, 2013


if isa(varargin{1},'GRIDobj')
    DEM = varargin{1};
    IX = SWATHobj2GRIDobj(SW,DEM,'ix');
    DX = SWATHobj2GRIDobj(SW,DEM,'distx');
    OUT = DEM; OUT.Z(:) = nan;
    varargin(1) = [];
    mkgrid = 1;
else
    OUT = SW; OUT.Z(:) = nan;
    mkgrid = 0;
end

p = inputParser;
p.FunctionName = 'baselevel';
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addParamValue(p,'smooth',5e2,@(x) isnumeric(x))
addParamValue(p,'values',5e2,@(x) isnumeric(x))
parse(p,SW,varargin{:});
smooth = p.Results.smooth;
values = p.Results.values;

if ~isempty(values)
    ix = find(isnan(values(:,1)));
    if ~isempty(ix) % multiple channel segments
        ix1 = ix(1:end)+1;
        ix2 = [ix(2:end)-1;length(values(:,1))];
    else
        error('Provided ''values'' have unknown format.')
    end
end

D = [];

for i = 1 : length(SW.xy0)
    
    d = SW.distx{i};
    d0 = SW.zd0{i}(:,2);
    
    if ~isempty(values)
        z0 = values(ix1(i):ix2(i),1);
        if size(values,2)==2
            d0 = values(ix1(i):ix2(i),2);
        end
    elseif (mkgrid)
        x0 = SW.xy0{i}(:,1);
        y0 = SW.xy0{i}(:,2);
        % If the SWATHobj is derived from a STREAMobj, then the following
        % conversion should result in no residuals
        [ix0,~] = coord2ind(DEM,x0,y0);
        z0 = DEM.Z(ix0);
    else
        z0 = SW.zd0{i}(:,1);
    end
    
    min_z = interp1(d0,z0,d);
%    min_z = min(SW.Z{i},[],1); % across-swath minimum
    
    % make sure to have a value for every point along the SWATHobj
    ix = ~isnan(min_z);
    zi = interp1(d(ix),min_z(ix),d);
    
    % there might still be NaNs at the ends that could blow up during smoothing
    if any(isnan(zi))
        ix = ~isnan(zi);
        zi = interp1(d(ix),zi(ix),d,'pchip','extrap');
    end
    
    % smooth the elevations
    if smooth>0
        try
            n = round(smooth/SW.dx);
            b = ones(n,1)./n;
            zf = filtfilt(b,1,double(zi));
        catch
            warning('Selected smoothing length too long for channel segment. Will use longest smoothing length possible.');
            n = round((length(zi)/3)-1);
            b = ones(n,1)./n;
            zf = filtfilt(b,1,double(zi));
        end
    else
        zf = zi;
    end
    
    % interpolat base-level grid with smooth elevations and distances
    if (mkgrid)
        ix_grid = find(IX.Z==i);
        OUT.Z(ix_grid) = interp1(d,zf,DX.Z(ix_grid));
    else
        OUT.Z{i} = repmat(zf,length(SW.disty{1}),1);
    end
    
    % append nan-separated list with d and z values
    D = [D; [nan,nan]; [d,zf]];
    
end


varargout{1} = OUT;
if nargout>1
    varargout{2} = D;
end

