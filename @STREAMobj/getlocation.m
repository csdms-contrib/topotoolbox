function varargout = getlocation(S,d0,varargin)

%GETLOCATION Get locations along a stream network
%
% Syntax
%     
%     [x,y,value]   = getlocation(S,val)
%     [x,y,value]   = getlocation(S,val,'value',S.distance)
%     [P,value]     = getlocation(...,'output','PPS','z',DEM)
%     [xc,yc,value] = getlocation(...,'output','xcyc')
%     [ix,value]    = getlocation(...,'output','ix')
%     MS            = getlocation(...,'output','mapstruct')
%     MP            = getlocation(...,'output','mappoint')
%     GP            = getlocation(...,'output','geopoint')
%     [x,y,value]   = getlocation(...,'output','cell')
%
% Description
%
%     getlocation returns the coordinates (or a data model that stores 
%     coordinates) where values measured along a stream network S have a
%     specified value val. For example, the command
%     [x,y,value]   = getlocation(S,val) or 
%                     getlocation(S,val,'value',S.distance)
%     returns the coordinates x and y at which the stream network has a
%     distance val from the outlet. Alternatively, the command
%     [x,y,value]   = getlocation(S,1000,'value',DEM) 
%     returns the coordinates where elevations along a stream network are
%     1000 m.
%     
% Input arguments
%
%     S    STREAMobj
%     val  scalar or vector of values (e.g. distance from the outlet or 
%          elevation), or anonymous function
%          'value' can also be a function such as @mean.
%          The command 
%          [x,y,value]   = getlocation(S,@mean,'value',DEM)
%          returns the locations where the stream network attains the mean 
%          elevation along the stream network.
%    
%     Parameter name/value pairs
%
%     'value'   node-attribute list or GRIDobj
%               By default, 'value' is S.distance, i.e. the distance from
%               the outlet. 
%     'output'  {'xy'},'xcyc','ix','PPS','mappoint','geopoint','mapstruct'
%               'xy'        exact coordinates (linearly interpolated)
%               'xcyc'      coordinates at cell centers
%               'ix'        linear indices of cells
%               'PPS'       instance of PPS
%               'mappoint'  mappoint
%               'geopoint'  geopoint (if S has a projected coordinate
%                           system)
%               'mapstruct' mapping structure array
%               'cell'      a cell array of points is returned for each
%                           value in val.
%     'z'       adds an elevation attribute to PPS (only applicable if
%               output is PPS)
%
% Output arguments
%
%     ...       see 'output'
%     value     same as val, but repeated for each location
%
% Example 1: Create map with points equally spaced from the outlet
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     [x,y,val] = getlocation(S,[0:5000:max(S.distance)]);
%     plot(S)
%     hold on
%     scatter(x,y,20,val/1000,'filled')
%     h = colorbar;
%     h.Label.String = 'Distance from outlet [km]';
%
% Example 2: Create a map with points where streams intersect with contour
%            lines
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     [xcon,ycon,zcon] = contour(DEM,10);
%     plot(xcon,ycon,'color',[.7 .7 .7])
%     P = getlocation(S,unique(zcon),'value',DEM,'output','PPS');
%     [x,y] = getlocation(S,unique(zcon),'value',DEM,'output','xy');
%     hold on
%     plot(P)
%     plot(x,y,'.','color',[.7 .7 .7])
%     hold off
%
% Example 3: Animation of upstream migrating knickpoints
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S,1);
%     A = flowacc(FD);
%     c = chitransform(S,A,'mn',0.4);
%     [xc,yc,val] = getlocation(S,linspace(min(c),max(c),500),...
%                               'value',c,'output','cell');
%     plot(S,'color',[.6 .6 .6])
%     axis image
%     hold on
%     h = plot(xc{1},yc{1},'ok','MarkerFaceColor',[.6 .6 .6]);
%     for r = 1:numel(val)
%         set(h,'XData',xc{r},'YData',yc{r});
%         drawnow
%     end
%     hold off
%         
%            
% See also: STREAMobj, STREAMobj/distance, STREAMobj/getvalue
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. November, 2020

% Valid output
validoutput = {'xy','xcyc','ix','PPS','mappoint',...
               'geopoint','mapstruct','cell'};

% Parse input
p = inputParser;
addParameter(p,'value',S.distance,@(x) isa(x,'GRIDobj') || isnal(S,x));
addParameter(p,'output','xy')
addParameter(p,'z',[])
parse(p,varargin{:});

% Validate output
output = validatestring(p.Results.output,validoutput);

% Handle value grid or node-attribute list
d = p.Results.value;
if isa(d,'GRIDobj')
    d = getnal(S,d); 
end

% Is val an function handle?
if isa(d0,'function_handle')
    d0 = d0(d);
end

% Outlets need special treatment
outlet = streampoi(S,'outlet','logical');

% Call getlocation_sub for all elements in d0
[x,y,c] = cellfun(@getlocation_sub, num2cell(d0),'UniformOutput',false);


switch lower(output)
    case 'cell'
    otherwise     
        x = vertcat(x{:});
        y = vertcat(y{:});
        c = vertcat(c{:});
end

% Prepare output
switch lower(output)
    case {'xy','xcyc'}
        varargout{1} = x;
        varargout{2} = y;
        varargout{3} = c;
    case 'ix'
        ix = coord2ind(GRIDobj(S,'logical'),x,y);
        
        varargout{1} = ix;
        varargout{2} = c;
    case 'mapstruct'
        varargout{1} = struct('Geometry','Point',...
            'X',num2cell(x),'Y',num2cell(y),'value',num2cell(c));
    case 'pps'
        if isempty(p.Results.z)
            varargout{1} = PPS(S,'PP',[x y]);
        else
            varargout{1} = PPS(S,'PP',[x y],'z',p.Results.z);
        end
        varargout{2} = c;
        
    case 'mappoint'
        varargout{1} = mappoint(x,y,'value',c);
    case 'geopoint'
        [lat,lon] = minvtran(S.georef.mstruct,x,y);
        varargout{1} = ygeopoint(lat,lon,'value',c);
    case 'cell'
        varargout{1} = x;
        varargout{2} = y;
        varargout{3} = cellfun(@(x) x(1),c,'UniformOutput',true);
end
%% --------------------------------------
function [x,y,c] = getlocation_sub(d0)
f = (d0-d(S.ixc)) ./ (d(S.ix)-d(S.ixc));
I = (f>0) & (f <=1) | (outlet(S.ixc) & (f == 0));

switch lower(output)
    case {'xcyc','ix','pps'} 
            f = round(f);
end
x = S.x(S.ixc(I)) + f(I).*(S.x(S.ix(I))-S.x(S.ixc(I)));
y = S.y(S.ixc(I)) + f(I).*(S.y(S.ix(I))-S.y(S.ixc(I)));

c = repmat(d0,size(x));
end
end
