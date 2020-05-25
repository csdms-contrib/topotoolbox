function DOUT = shrink(DIN,FD,varargin)
%SHRINK   Shrink divide network
%
% Syntax
%
%     DOUT = shrink(DIN,FD)
%     DOUT = shrink(DIN,FD,distance)
%     DOUT = shrink(DIN,FD,type)
%     DOUT = shrink(DIN,FD,type,val)
%
% Description
%
%     SHRINK reduces the extent of a divide network based on
%     either entire segments, the divide order, or distance.
%     Note that shrinking based on segments and divide order
%     does not change the length of segments, whereas shrinking
%     based on distance will likely do so. Input divide objects 
%     need to be sorted. 
%
% Input
%
%     DIN       instance of class DIVIDEobj
%     distance  threshold distance used for shrinking divide networks
%     type      character array indicating shrinking type
%               ('segment' (default), 'order', 'distance')
%     val       threshold value used for shrinking. Corresponds
%               to the number of segments to remove, the highest
%               order, or the maximum distance
%
% Output
%
%     DOUT      instance of class DIVIDEobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,ST);
%     D2 = shrink(D,FD,1000);
%     subplot(2,1,1)
%     plot(D,'color','k')
%     axis image
%     subplot(2,1,2)
%     plot(D2,'color','r')
%     axis image
%
% See also: DIVIDEobj, DIVIDEobj/plot
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020


if not(DIN.issorted)
    error('Sorted divides needed. Use function "sort" first.')
end

if nargin==2
    type = 'segment';
elseif nargin==3 && ischar(varargin{1})
    type = varargin{1};
    validatestring(type,{'segment','order'});
    val = 0;
elseif nargin==3 && isscalar(varargin{1})
    type = 'distance';
    val = varargin{1};
elseif nargin==4 && isscalar(varargin{2})
    type = varargin{1};
    validatestring(type,{'distance','segment','order'});
    val = varargin{2};
else
    error('Unknown call to function. See help.')
end

if ~strcmp(DIN.ordertype,'topo') && strcmp(type,'segment')
    warning(['Ordertype is "%s". Shrinking by segment',...
        'may produce unexpected results'],type);
end


DOUT = DIN;


switch type
    case 'segment'
        if val==0
            val = 1;
        end
        if not(strcmp(DOUT.ordertype,'topo'))
            tDOUT = divorder(DOUT,'topo');
        else
            tDOUT = DOUT;
        end
        dix = tDOUT.order>val | isnan(tDOUT.IX);
        
    case 'order'
        if val==0
            val = min(unique(DIN.order));
        end
        dix = DIN.order>val | isnan(DIN.IX);
        
    case 'distance'
        if isempty(DIN.distance)
            DIN = divdist(DIN);
        end
        dix = DIN.distance>val | isnan(DIN.IX);
        
end
fix = find(not(dix))+1;
inan = isnan(DOUT.IX(fix));
dix(fix(inan)) = false;

IX = DOUT.IX(dix);
M = onl2struct(IX);
nr = [M.nr];
ix = nr>1;
DOUT.IX = vertcat(M(ix).IX);
DOUT = divnet(DOUT,FD);
DOUT = sort(DOUT);

if ~isempty(DIN.distance)
	DOUT = divdist(DOUT);
end
if ~isempty(DIN.order)
	DOUT = divorder(DOUT,DIN.ordertype);
end 

end
