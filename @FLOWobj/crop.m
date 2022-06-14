function [FD,W] = crop(FD,mask)

%CROP crop an instance of FLOWobj
%
% Syntax
%
%     FDc = crop(FD,MASK)
%     [FDc,W] = crop(FD,MASK);
%
% Description
%
%     This function crops an instance of FLOWobj to the extent indicated by
%     logical values in MASK. MASK is returned by the function GRIDobj/crop
%     as second output argument. 
%
% Input arguments
%
%     FD    FLOWobj
%     MASK  logical GRIDobj as returned by GRIDobj/crop. FD and MASK must
%           spatially align (see validatealignment)
%
% Output arguments
%
%     FDc   cropped FLOWobj
%     W     GRIDobj for a weighted flow accumulation that preserves the 
%           incoming fluxes.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     [DEMc,MASK] = crop(DEM,DEM>1300,nan);
%     FDc = crop(FD,MASK);
%     A   = flowacc(FDc);
%     imageschs(DEMc,log(A))
%
% 
% See also: GRIDobj/crop, GRIDobj/pad
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 2. February, 2021

validatealignment(FD,mask)

if ~isa(mask,'GRIDobj')
    temp = GRIDobj(FD,'logical');
    temp.Z = mask;
    mask = temp;
end

if nargout == 2
    A = flowacc(FD);
    W = GRIDobj(FD,'single') + 1;
    % Get those edges directing inwards into the mask
    I     = (~mask.Z(FD.ix)) & (mask.Z(FD.ixc));
    for r = 1:numel(FD.ix)
        if I(r)
            W.Z(FD.ixc(r)) = W.Z(FD.ixc(r))+A.Z(FD.ix(r));
        end
    end
   
    W = crop(W,mask);
end


I = all([mask.Z(FD.ix) mask.Z(FD.ixc)],2);
FD.ix = FD.ix(I);
FD.ixc = FD.ixc(I);
if ismulti(FD)
    FD.fraction = FD.fraction(I);
end

% map old linear indices to new linear indices
[I,mask] = crop(GRIDobj(FD,'logical'),mask);
[r,c] = find(mask);
mask.Z(min(r):max(r),min(c):max(c)) = true;
mask.Z   = uint32(mask.Z);
mask.Z(mask.Z>0) = uint32(1:nnz(mask.Z));
FD.ix  = mask.Z(FD.ix);
FD.ixc = mask.Z(FD.ixc);

FD.size = I.size;
FD.refmat = I.refmat;
FD.georef = I.georef;
FD.fastindexing = false;

