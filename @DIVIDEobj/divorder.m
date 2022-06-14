function DOUT = divorder(DIN,varargin)
%DIVORDER   Assign order to divide segments
%
% Syntax
%
%     D2 = divorder(D)
%     D2 = divorder(D,type)
%
% Description
%
%     DIVORDER assigns to each segment in the divide network a divide
%     order, analogous to stream ordering. The available ordering schemes
%     are 'strahler','shreve', and 'topo'. In the Strahler ordering scheme,
%     the order increases by one if the joining divide segments have the 
%     same order, otherwise it remains at their maximum order. In the 
%     Shreve ordering scheme, the resulting divide order is the sum of 
%     those of the joining divide segments, and in the Topo ordering 
%     scheme, divide orders increase by one at each junction.
%
% Input
%
%     D         instance of class DIVIDEobj
%     type      character array indicating the ordering scheme
%               ('strahler', 'shreve', or 'topo' (default))
%
% Output
%
%     D2        instance of class DIVIDEobj
%
% Example
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,ST);
%     D = sort(D);
%     D = divorder(D,'strahler');
%     D2 = divorder(D,'topo');
%     subplot(2,1,1)
%     plot(D,'color','k')
%     title('Strahler')
%     axis image, box on
%     subplot(2,1,2)
%     plot(D2,'color','k')
%     title('Topo')
%     axis image, box on
%
%
% See also: DIVIDEobj, DIVIDEobj/sort
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020


if not(DIN.issorted)
    error('Divides need to be sorted. Use function SORT first.')
end

DOUT = DIN;

if nargin>1
    type = lower(varargin{1});
else
    type = 'topo';
end

% divide order function
switch type
    case 'strahler'
        dofct = @(do)max([min(do)+1;do]);
    case 'shreve'
        dofct = @(do)sum(do);
    case 'topo'
        dofct = @(do)max(do)+1;
    otherwise
        error('Unknown type')
end

% Prepare
M = onl2struct(DIN.IX);
[M.do] = deal([]);

st = vertcat(M.st);
st1 = st(:,1);
st2 = st(:,2);
[~,lib] = ismember(st2,st1);

% set divide order (1) for endpoints
ixep = ismember(st1,DIN.ep);
[M(ixep).do] = deal(1);
[M.checked] = deal(false);
ixjct = not(ixep);

% loop over divide segments
for i = 1 : length(M)
    ixc = lib(i);
    if ixjct(i)
        M(i).do = dofct(M(i).do).*ones(size(M(i).IX));
        %M(i).do(end) = NaN;
        M(i).checked = true;
    else
        M(i).do = ones(size(M(i).IX));
        %M(i).do(end) = NaN;
        M(i).checked = true;
    end
    if ixc>0
        if not(M(ixc).checked)
            M(ixc).do = [M(ixc).do; max(M(i).do)];
        end
    end
end

DOUT.order = vertcat(M.do);
DOUT.ordertype = type;

end
