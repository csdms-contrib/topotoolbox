function direc = tribdir(S)

%TRIBDIR direction of inflow of tributary
%
% Syntax
%
%     direc = tribdir(S)
%
% Description
%
%     tribdir computes whether a stream segment confluences from the right
%     or left direction. 
%
% Input arguments
%
%     S      STREAMobj
%     
% Output arguments
%
%     direc  node-attribute list. Stream segments have values of 1, if they
%            enter confluences from the right, and -1 if they come from the
%            left. The direction of stream segments with a value of 0 
%            cannot be determined (e.g. single rivers connected to the edge
%            of the DEM).
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = klargestconncomps(S);
%     direc = tribdir(S);
%     subplot(1,2,1)
%     plotc(S,direc)
%     axis image
%     box on
%     subplot(1,2,2)
%     plotdz(S,DEM,'color',direc)
%     
% See also: STREAMobj/orientation, STREAMobj/curvature
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 24. May, 2019


I = streampoi(S,'confl','logical');

II = I(S.ixc);

IX = [S.ix(II) S.ixc(II)];

[II,loc] = ismember(IX(:,2),S.ix);
IX(:,3) = zeros(size(IX,1),1);
IX(II,3) = S.ixc(loc(II));

% confluences at the outlet are somewhat special
oc = ~II;    
n  = numel(S.IXgrid);
c  = accumarray(IX(oc,2),IX(oc,1),[n 1],@sum);
IX(oc,3) = c(IX(oc,2))-IX(oc,1);

% there are triple confluences at DEM egdes. These are exceptions that are
% deleted. Yuc...
delet = sum(sparse(IX(:,2),IX(:,1),true,n,n),2) > 2;
delet = delet(IX(:,2)) & oc;
IX(delet,:) = [];
oc(delet,:) = [];


d = (S.y(IX(:,1))-S.y(IX(:,2))).*(S.x(IX(:,3))-S.x(IX(:,1))) - ...
        (S.x(IX(:,1))-S.x(IX(:,2))).*(S.y(IX(:,3))-S.y(IX(:,1)));
    

dd = accumarray(IX(:,2),d,[n 1],@max);
md = dd(IX(:,2));
I  = d < md;
I  = +I;
% I(I==0) = -1; 

dd = accumarray(IX(:,2),d,[n 1],@min);
md = dd(IX(:,2));
II = d > md;
I(II) = -1;

I(oc) = -I(oc);

direc = zeros(n,1);
direc(IX(:,1)) = I;

for r = numel(S.ix):-1:1
    if direc(S.ix(r))==0
        direc(S.ix(r)) = direc(S.ixc(r));
    end
end
direc(isnan(direc)) = 0;        
