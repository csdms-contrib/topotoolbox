function DOUT = removeshortdivides(DIN,FD,d)
%REMOVESHORTDIVIDES Remove short first-order divides
%
% Syntax
%
%     DOUT = removeshortdivides(DIN,FD,d)
%
% Description
%
%     Divide networks sometimes include short first order divides
%     that one may want to exclude from further analysis.
%     REMOVESHORTDIVIDES enables to remove first-order divides 
%     with a length equal or less to d.
%
% Input arguments
%
%     DIN   divides (class DIVIDEobj)
%     d     length of first-order divides to be removed in map
%           units
%
% Output arguments
%
%     DOUT  divides (class DIVIDEobj)
%
% Example
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     ST = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,ST);
%     D = divorder(D,'topo');
%     plot(D,'color','r','limits',[0 inf])
%     hold on
%     D = removeshortdivides(D,FD,500);
%     plot(D,'color','k','limits',[0 inf])
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020

if not(DIN.issorted)
    error('Sorted divides needed. Use function "sort" first.')
end

DOUT = DIN;
tDOUT = divorder(DIN,'topo');

% get indices of separating NaNs
ix = find(isnan(DIN.IX));
% get length of the following divides
dix = diff(ix);
% get divides with length less then d
fix = find(dix<d./DIN.cellsize);
ND = DIN.IX;
for k = 1 : length(fix)
    ND(ix(fix(k))+1:ix(fix(k))+dix(fix(k))) = -1;
end
II = not(ND<0) | tDOUT.order>1;
DOUT.IX = DIN.IX(II);
DOUT = divnet(DOUT,FD);
DOUT = sort(DOUT);

if ~isempty(DIN.distance)
	DOUT = divdist(DOUT);
end
if ~isempty(DIN.order)
	DOUT = divorder(DOUT,DIN.ordertype);
end 

end
