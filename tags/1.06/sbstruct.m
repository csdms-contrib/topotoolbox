function  [S,BASINID] = sbstruct(M,ix,siz)

% create structure array for sub-basin analysis
%
% Syntax
%
%     [S,SB] = sbstruct(M,ix,siz)
%
% Description
%
%     sbstruct creates a structure that can subsequently be used for
%     sub-basin analysis (sbprops). 
%
% Input
%
%     M     single flowdirection matrix
%     ix    vector with linear indices of subbasin gauges
%     siz   size of the DEM
%
% Output
%
%     S     structure array (use as input for sbplot and sbprops)
%           
%           .SBIx: gauge and subbasin id
%           .OutletIdx: linear index auf gauges (=ix)
%           .NextIx: id of next downstream subbasin gauge (if zero then 
%           there is no downstream gauge)
%           .SBIdxList: Linear index of pixels draining into closest
%           subbasin gauge
%           .PreviousIx: all upslope gauges
%           .Adj: Adjacency matrix of sub-basins (rows: out, cols: in)
%     SB    label matrix containing subbasins
%
% Example 
%
%     (see sbplot)
%
%
% See also: COORD2IND, SBPLOT, SBPROPS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 30. October, 2009





error(nargchk(2, 3, nargin))

if nargout == 2 && nargin == 2
    error('TopoToolbox:incorrectinput',...
        'sbstruct requires siz as third input argument \nwhen called with two output arguments.')
end


M = spones(M);
% do we have a multiple flow direction matrix?
if any(sum(M,2)>1);
    error('TopoToolbox:incorrectinput',...
        'single flow direction must be supplied')
end

nrc = size(M,1);

if nargout == 2 && nrc~=prod(siz);
    error('TopoToolbox:incorrectinput',...
        'prod(siz) must equal to size(M,1) and size(M,2)')
end

if any(ix>nrc)
    error('TopoToolbox:incorrectinput',...
        'ix exceeds matrix dimensions')  
end

if numel(ix) ~= numel(unique(ix));
    error('elements in ix must be unique')
end

% force column vector
ix = ix(:);

% nr of outlets
nroutlets   = numel(ix);
% vector with station index increments
S.SBIx      = (1:nroutlets)';
% outlet indices
S.OutletIdx = ix;

% nr of cells in the dem
nrc         = size(M,1);  

% get connectivity between sub-basin outlets
% find next subbasin outlet
[ixun,ixix] = find(M(:,ix));
OUTLETS     = zeros(nrc,1);
OUTLETS(ix) = S.SBIx;
Mtemp       = sparse(ix(ixix),ixun,1,nrc,nrc)';
NextIx     = (speye(nrc)-(M>Mtemp))\(Mtemp*OUTLETS);
S.NextIx   = NextIx(ix);


% find linear indices of pixels belonging to next subbasins
[ixix,ixdn] = find(M(ix,:));
BASINID     = (speye(nrc)-(M>sparse(ix(ixix),ixdn,1,nrc,nrc)))\OUTLETS;
I           = BASINID > 0;
S.SBIdxList = accumarray(BASINID(I),find(I),[nroutlets 1],@(x) {x});

% find upslope basins 
I = S.NextIx ~= 0;
S.Adj = sparse(S.SBIx(I),S.NextIx(I),1,nroutlets,nroutlets);
IN = zeros(nroutlets,1);
for outlet = 1:nroutlets;
    IN(S.SBIx(outlet)) = 1;
    S.PreviousIx{outlet,1} = find(((speye(nroutlets)-S.Adj)\IN)>IN);
    IN(:) = 0;
end

% create array with BASINIDs if requested as second output argument
if nargout==2
    BASINID = reshape(BASINID,siz);
end


