function [ic,icd] = routegeodesic(DEM,type)

% use least cost path to route through flats of a digital elevation model
%
% Syntax
%
%     [IXf,IXn] = routegeodesic(dem,type)
%
% Description
%
%     routegeodesic is a subroutine used by some of the flowdirection
%     algorithms in the toolbox. routegeodesic routes through flats in a
%     very natural way by going along the centerline of flat areas. The
%     user may choose between single and multiple flow routing. The output
%     are vectors containing linear indices of connected cells along the
%     flow direction. IXn are the downstream neighbors of IXf.
%
% Input
% 
%     dem       digital elevation model
%     type      string for flow routing method ('single' (default) or
%               'multi')
%
% Output
%     
%     IXf       vector containing linear indices of cells in flats
%     IXn       vector containing linear indices of downstream neighbors
%               of IXf
%
%
% See also: CROSSFLATS, ROUTEFLATS
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 7. October, 2011


dem = DEM.Z;

siz = size(dem);

[I,SILLS] = identifyflats(DEM);
I = I.Z;
SILLS = SILLS.Z;

if any(I(:))
    
    % establish the connectivity between sills and flats
    [SillPixel,PreSillPixel] = ixneighbors(SILLS,SILLS);
    II = (dem(SillPixel) == dem(PreSillPixel)) & I(PreSillPixel) & SILLS(SillPixel);
    SillPixel    = SillPixel(II);
    PreSillPixel = PreSillPixel(II);
    % use PreSillPixel later as marker for distance transform
    
    % calculate weight matrix for geodesic distance transform
    D  = bwdist(~I);
    S  = regionprops(I,D,'MaxIntensity','PixelIdxList');
    for r = 1:numel(S);
        D((S(r).PixelIdxList)) = S(r).MaxIntensity - D((S(r).PixelIdxList)) + 1;
    end
    D(~I)    = inf;
    
    C  = graydist(D,PreSillPixel,'quasi-euclidean');
    
    ic = find(I);
    
    neighs = [-1 -1+siz(1) siz(1) siz(1)+1 1 1-siz(1) -siz(1) -1-siz(1)];
    distan = repmat([1 sqrt(2)],1,4);
    
    switch type
        case 'single'
            g = zeros(numel(ic),1);
            icd = zeros(numel(ic),1);
            for r = 1:8;
                ict = ic+neighs(r);
                gt  = (C(ic)-C(ict))/distan(r);
                UD  = gt>g;
                icd(UD) = ict(UD);
                g(UD)   = gt(UD);
            end
            II = icd == 0;
            icd(II) = [];
            ic(II) = [];
        case 'multi'
            icd = [];
            ic2 = [];
            for r = 1:8;
                DSN = ic+neighs(r);
                II  = C(ic)-C(DSN) > 0;
                icd = [icd; DSN(II)];
                ic2 = [ic2; ic(II)];
                
            end
            ic = ic2;
    end
    
    
    ic  = [ic;  PreSillPixel];
    icd = [icd; SillPixel];
else
    ic = [];
    icd = [];
end
