function [IXf,IXn] = routeflats(dem,type)

% route through flats of a digital elevation model
%
% Syntax
%
%     [IXf,IXn] = routeflats(dem,type)
%
% Description
%
%     routeflats is a subroutine used by some of the flowdirection 
%     algorithms in the toolbox. routeflats recursively creates flow paths
%     through flat terrain. The user may choose between single and multiple
%     flow routing. The output are vectors containing linear indices of
%     connected cells along the flow direction. IXn are the downstream
%     neighbors of IXf.
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
% See also: CROSSFLATS
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 17. June, 2010



% type: 'single' or 'multi'
type     = lower(type);

% floating-point relative accuracy
fpra     = eps(dem) * 2;
fillval  = min(dem(:))-max(fpra(:));

% kernel (8-Neighborhood)
nhood    = ones(3,3);
siz      = size(dem);

% handle NaNs
log_nans = isnan(dem);
if any(log_nans(:));
    flag_nans = true;
else
    flag_nans = false;
end

% if there are nans.
if flag_nans    
    dem(log_nans) = -inf;
end

% identify flats
% flats: logical matrix with true where cells don't have lower neighbors
if flag_nans
    flats = imerode(dem,nhood) == dem & ~log_nans;
else
    flats = imerode(dem,nhood) == dem;
end

% remove regional minima
flats(imregionalmin(dem)) = false;

% if there are nans.
if flag_nans    
    dem(log_nans) = fillval;
end

% remove flats at the border
flats(:,[1 2 end-1 end])  = false;
flats([1 2 end-1 end],:)  = false;

% any flats there? If yes, start recursive flat routing
if any(flats(:));
    II = true;
    
    switch type
        case 'single'
            [IXf,IXn] = routeflatssingle(flats);
        
        case 'multi'    
            [IXf,IXn] = routeflatsmulti(flats);
            
        otherwise
            IXf = [];
            IXn = [];
    end
    
    if flag_nans
        IXf(log_nans(IXn)) = [];
        IXn(log_nans(IXn)) = [];
    end
else
    IXf = [];
    IXn = [];
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
function [ic,icd] = routeflatssingle(flats)
% recursive routing through flats using single flow direction
%
% [ic,icd] = routeflatssingle(flats)


% linear indices of flat cells
ic   = find(flats);
% number of flat cells
nrcf = numel(ic);
% create nan array that will be filled with downstream
% neigbors of indices in ic in
% the subsequent lines of codes 
icd = nan(nrcf,1);

% inner border in flats
% finline   = bwmorph(flats,'remove',8); 
finline   = bwperim(flats,8); 

% first create connectivity to downstream cells 
Igive   = finline(ic);
IXgive  = find(Igive);

% find neighbors of 
IXn = neighs(ic(Igive));
I2  = ~ismember(IXn,ic(Igive)) & ~flats(IXn) & bsxfun(@le,dem(IXn),dem(ic(Igive)));
% cells to which an flow direction can be assigned
% do have at least one downstream neighbor (I2)
Igive(IXgive(~any(I2,2))) = false;

I2  = I2';
I2  = I2 & cumsum(I2,1) == 1;

IXn = IXn';

icd(Igive) = IXn(I2(:));

% indicate the cells as nonflat
flats(ic(Igive)) = false;

% another logical array that indicates if a flat cell receives from
% upstream neighbor
Ireceive = false(nrcf,1);

% iteratively remove non-giving cells
while any(~Igive)
    
    % find cells that give but not receive
    Itemp = Igive & ~Ireceive;
    IXn   = ic(Itemp);
   
    % find neighbors of these cells
    IXf = neighs(IXn);
    % which of these neighbors are flat cells
    I2  = flats(IXf);
    % which IXn do have flat neighbors at all
    
    II  = any(I2,2);
    
    IXn = IXn(II);
    IXf = IXf(II,:);
    I2  = I2(II,:);
    
    IXn = repmat(IXn,1,8)';
    
    IXf = IXf';
    I2  = I2';
    
    IXn = IXn(I2(:));
    IXf = IXf(I2(:));
    
    [IXf,m] = unique(IXf);
    IXn     = IXn(m); 
    
    I  = ismembc(ic,IXf);
    Igive(I) = true;
    icd(I) = IXn;
    Ireceive(Itemp) = true;
    
    flats(IXf) = false;   
    
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ic,icd] = routeflatsmulti(flats)
% recursive routing through flats using multiple flow direction
%
% [ic,icd] = routeflatsmulti(flats)
%
%
%
nrflats = sum(flats(:));

ic = cell(nrflats,1);
icd = cell(nrflats,1);

counter = 1;
% recursive routing through flats
while II
    
    % inner border in flats
    finline   = bwmorph(flats,'remove',8);
    % outline around flats
    outline   = imdilate(flats,nhood) & ~flats;
    
    % outline
    IXn = find(outline);
    IXf = neighs(IXn);
    I   = finline(IXf) & bsxfun(@le,dem(IXn),dem(IXf));
    sI  = sum(I,2);
    I2  = sI>0 & sI<8;
    
    % inner border cells that don't have neighbors of
    % equal or lower elevation are sinks
    flats(IXf(sI==8,:)) = false;
    
    IXf = IXf(I2,:);
    IXn = IXn(I2);
    I   = I(I2,:);
    
    IXn  = spdiags(IXn,0,numel(IXn),numel(IXn))*I;
    IXn  = IXn(I);
    IXf  = IXf(I);
    
    flats(IXf(:)) = false;
    
    ic{counter} = full(IXf(:));
    icd{counter} = full(IXn(:));
    
    if any(flats(:))
        counter = counter+1;
    else
        II = false;
        ic = cell2mat(ic);
        icd = cell2mat(icd);
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IXn = neighs(IXf)
% get neighbors of a cell in a matrix

    ix = randperm(8);
    IXn  = zeros(size(IXf,1),8);
    IXn(:,ix(1)) = IXf+1;
    IXn(:,ix(2)) = IXf-1;
    IXn(:,ix(3)) = IXf+siz(1);
    IXn(:,ix(4)) = IXf-siz(1);

    IXn(:,ix(5)) = IXf+1+siz(1);
    IXn(:,ix(6)) = IXf-1+siz(1);
    IXn(:,ix(7)) = IXf-1-siz(1);
    IXn(:,ix(8)) = IXf+1-siz(1);
end


end
