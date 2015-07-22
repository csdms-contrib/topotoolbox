function z = mincosthydrocon(S,DEM,method,fillp)

% minimum cost hydrological conditioning
%
% Syntax
%
%     z = mincosthydrocon(S,DEM)
%     z = mincosthydrocon(S,DEM,method)
%     z = mincosthydrocon(S,DEM,'interp',fillp)
%
% Description
%
%     mincosthydrocon implements different methods to correct channel
%     length profiles prone to data artifacts. These artifacts often lead
%     to the unwanted behavior that profiles are not monotically decreasing
%     from the source to the outlet.
%     
%     The method 'minmax' locally chooses between carving and filling
%     thereby minimizing the sum of modifications (sum(z-z_carve) +
%     sum(z_fill-z)) required to obtain a hydrologically correct channel
%     length profile.
%     
%     The method 'interp' linearly interpolates in longitudinal channel
%     sections that are hydrologically incorrect, i.e., sections that must
%     either be carved or filled.
%
% Input arguments
%
%     S         STREAMobj
%     DEM       digital elevation model (GRIDobj)
%     method    'minmax' (default) or 'interp'
%     fillp     scalar between 0 and 1, determines the degree of "filling".
%               Values close to 0 lead to a higher degree of carving.
%
% Output arguments
%
%     z         node attribute list (nal) that contains the modified 
%               elevation values along the stream network S.
%
% Example
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S  = STREAMobj(FD,flowacc(FD)>500);
%     z = mincosthydrocon(S,DEM);
%     plotdz(S,DEM,'color','k');
%     hold on
%     plotdz(S,z,'color','r');
%     z = mincosthydrocon(S,DEM,'interp');
%     plotdz(S,z,'color','b');
%
% See also: STREAMobj, STREAMobj/plotdz
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 24. June, 2014


narginchk(2,4);

if nargin == 2;
    method = 'minmax';
elseif nargin == 3;
    method = validatestring(method,{'minmax', 'interp'},'mincosthydrocon','method',3);
    fillp  = .5;
else
    method = validatestring(method,{'minmax', 'interp'},'mincosthydrocon','method',3);
    validateattributes(fillp,{'numeric'},{'scalar','>',0,'<=',1});
end

d      = S.distance;
% node attribute list (nal) with filled z
z      = DEM.Z(S.IXgrid);
z_fill = DEM.Z(S.IXgrid);

for r = numel(S.ixc):-1:1;
    z_fill(S.ix(r)) = max(z_fill(S.ix(r)),z_fill(S.ixc(r)));
end

% nal with carved z
z_carve = DEM.Z(S.IXgrid);
for r = 1:numel(S.ixc);
    z_carve(S.ixc(r)) = min(z_carve(S.ix(r)),z_carve(S.ixc(r)));
end

ICON = (z_carve-z)~= 0 | (z_fill-z)~= 0;

ix = S.ix;
ixc = S.ixc;

IEDGE = ICON(ix);% | ICON(ixc);
ix = ix(IEDGE);
ixc = ixc(IEDGE);

switch method
    case 'minmax'
        c_fill = z_fill-z;
        for r = 1:numel(ixc);
            c_fill(ixc(r)) = c_fill(ix(r)) + c_fill(ixc(r));
        end
        
        for r = numel(ixc):-1:1;
            c_fill(ix(r)) = max(c_fill(ix(r)),c_fill(ixc(r)));
        end
        
        c_carve = z-z_carve;
        for r = 1:numel(ixc);
            c_carve(ixc(r)) = c_carve(ix(r)) + c_carve(ixc(r));
        end
        
        for r = numel(ixc):-1:1;
            c_carve(ix(r)) = max(c_carve(ix(r)),c_carve(ixc(r)));
        end
        
        I = c_fill<c_carve;
        z(I) = z_fill(I);
        z(~I) = z_carve(~I);
    case 'interp'
%         for r = 1:numel(ix);
%             z_fill(ixc(r)) = max(z_fill(ix(r)),z_fill(ixc(r)));
%         end
%         for r = numel(ix):-1:1;
%             z_carve(ix(r)) = min(z_carve(ixc(r)),z_carve(ix(r)));
%         end
        
        d = zeros(size(z_fill));
        for r = numel(ix):-1:1;
            d(ix(r)) = d(ixc(r))+...
                sqrt((S.x(ixc(r))-S.x(ix(r)))^2 + (S.y(ixc(r))-S.y(ix(r)))^2);
        end
        d2 = d;
        for r = 1:numel(ix);
            d2(ixc(r)) = max(d2(ix(r)),d2(ixc(r)));
        end
        d = d./d2;
        
        z(ICON) = z_carve(ICON) + (z_fill(ICON) - z_carve(ICON)).*d(ICON)*fillp;
        
               
end


        

