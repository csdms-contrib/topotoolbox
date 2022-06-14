function DOUT = cleanedges(DIN,FD)
%DIVIDEobj   Delete divides on the edges of DEM
%
% Syntax
%
%     D2 = cleanedges(D,FD)
%
% Description
%
%     Divide objects that were derived from a DEM often contain
%     unreasonable divides that coincide with the edges of the input grid.
%     CLEANEDGES removes these divide segments and reorders the divide
%     network.
%
% Input arguments
%
%     D      instance of divide object (DIVIDEobj)
%     FD     instance of flow direction object (FLOWobj)
%
% Output arguments
%
%     D2      instance of DIVIDEobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,flowacc(FD)>1000);
%     D = DIVIDEobj(FD,S);
%     subplot(1,2,1)
%     plot(D)
%     D2 = cleanedges(D,FD);
%     subplot(1,2,2)
%     plot(D2)
%  
% See also: DIVIDEobj, DIVIDEobj/divnet
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: February 2019


nr = DIN.size(1);
nc = DIN.size(2);
r = [ones(nc,1);reshape([2:nr-1;2:nr-1],2*nr-4,1);ones(nc,1).*nr];
c = [(1:nc)';repmat([1;nc],nr-2,1);(1:nc)'];
iedge = sub2ind(DIN.size,r,c);

M = onl2struct(DIN.IX);

cleansegments = struct;
ct = 0;
for i = 1 : length(M)
    ix = M(i).IX;
    %thisdo = segments(i).do;
    ia = ismember(ix,iedge);
    if sum(ia)>0
        ix(ia) = NaN;
        % Remove redundant NaN
        I = isnan(ix);
        I = I(:);
        J = [I(2:end);true];
        IJ = logical(min([I J],[],2));
        ix = [ix(not(IJ),:);NaN];
        if isnan(ix(1))
            ix = ix(2:end);
        end
        %iy = [0;find(isnan(ix),1)];
        iy = [0;find(isnan(ix))];
        for k = 1 : length(iy)-1
            ct = ct+1;
            cleansegments(ct).ix = ix(iy(k)+1:iy(k+1));
        end
    else
        ct = ct+1;
        cleansegments(ct).ix = ix;
    end
end

DOUT = DIN;
DOUT.IX = vertcat(cleansegments.ix);
DOUT = divnet(DOUT,FD);
if DIN.issorted
    DOUT = sort(DOUT);
end
if not(isempty(DIN.distance))
    DOUT = divdist(DOUT);
end
if not(isempty(DIN.ordertype))
    DOUT = divorder(DOUT,DIN.ordertype);
end

end
