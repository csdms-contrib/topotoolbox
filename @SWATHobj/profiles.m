function [d,z,varargout] = profiles(SW,varargin)
% PROFILES obtains profiles from a SWATHobj at distinct positions
%
% Syntax
%
%    [d,z] = profiles(SW)
%    [d,z] = profiles(SW,'pn','pv',...)
%    [d,z,x,y] = profiles(SW,...)
%
%
% Description
%
%     Profiles is used to samples profiles from a SWATHobj in either the
%     along- or across-track direction
%
%
% Input arguments
%
%     SW     Swath profile object (Class: SWATHobj)
%
%     Parameter name value pairs:
%
%       'dist'      {'x'},'y'
%           specifies the direction along which the profiles run. 'x' 
%           means across-track and 'y' means along-track of the swath.
%
%       'step'      scalar {width of SWATHobj}
%           specifies the step between the profiles in map units.
%       
%
% Output
%
%     d     distance values
%     z     z-values 
%     x     x-coorindates
%     y     y-coordinates
%
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,A>100);
%     S = trunk(klargestconncomps(S,1));
%     [x,y] = STREAMobj2XY(S);
%     ix = ~isnan(x);
%     SW = SWATHobj(DEM,x(ix),y(ix),'width',3e3,'smooth',200);
%     [d,z,x,y] = profiles(SW,'dist','x','step',5e3);
%     figure,plot(SW), hold on, plot(x,y,'k--'), hold off
%     figure,plot(d,z)
%
%
% See also: SWATHobj
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: March, 2016




% Parse inputs
p = inputParser;
p.FunctionName = 'profiles';
addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
addParamValue(p,'dist','x',@(x) validatestring(x,{'x','y'}))
addParamValue(p,'step',SW.width,@(x) isnumeric(x))
parse(p,SW,varargin{:});

Z = SW.Z;

switch p.Results.dist
    case 'x'
        dx = SW.dx;
        dist_along = SW.distx;
        dist_across = SW.disty;
    case 'y'
        dx = SW.dy;
        dist_along = SW.disty;
        dist_across = SW.distx;
        Z = Z';
end

dix = round(p.Results.step/dx);
ix = 1:dix:length(dist_along);

xdata = cell(length(ix),1);
ydata = xdata;
for i = 1 : length(ix)
    v = Z(:,ix(i));
    ydata{i} = [v(:);nan];
    xdata{i} = [dist_across(:);nan];
end

if nargout>2
    xcoord = xdata;
    ycoord = xdata;
    for i = 1 : length(ix)
        xval = SW.X(:,ix(i));
        yval = SW.Y(:,ix(i));
        xcoord{i} = [xval(:);nan];
        ycoord{i} = [yval(:);nan];
    end
    varargout{1} = cell2mat(xcoord);
    varargout{2} = cell2mat(ycoord);
end

d = cell2mat(xdata);
z = cell2mat(ydata);



end

