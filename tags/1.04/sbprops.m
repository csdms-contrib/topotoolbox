function  P = sbprops(S,X,Y,agg,A)

% measure properties of sub-basins
%
% Syntax
%
%     P = sbprops(S,X,Y,agg,A)
%
% Description
%
%     sbprops(S,X,Y,agg,A) measures a set of properties for each sub-basin
%     in S for the matrix A. A must have same extent and size ans X and Y. 
%     The user can choose between two modes. agg = true or agg = false
%     (default). If agg == false, sbprops returns the properties of the 
%     sub-basins only. If agg == true, sbprops aggregates the properties
%     according to their connectivity.
%
% Input
%
%     S     sub-basin structure array as returned by sbstruct
%     X,Y   coordinate matrices of the DEM
%     agg   aggregation mode (false (default) or true)
%     A     Matrix with same size as X and Y, whose properties are to be
%           returned
%
% Output
%
%     P     output structure array. Currently sbprops supports following
%           parameters:
%           .Area
%           .Mean
%           .Sum
%           .Min'
%           .Max'
%           .Range' 
%
% Example
%
%     load exampleDEM
%     M = flowdir(X,Y,dem,'type','single');
%     IX = [12633 10683 8508 6227 8219 17263 17372]';
%     S = sbstruct(M,IX);
%     P = sbprops(S,X,Y,true,dem);
%     figure
%     scatter([P.Area],[P.Mean])
%     xlabel('Contributing area [m^2]')
%     ylabel('mean elevation in contributing area [m]')
%
%
%
% See also: SBPLOT, SBPROPS, REGIONPROPS, DRAINAGEDENSITY
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 12. Octobre, 2009



% cellsize
cs    = abs(Y(1)-Y(2));
% nr of subbasins
nSBIx = numel(S.SBIx);

% siz = size(X);

if nargin == 5
    pixvalmeas = true;
else
    pixvalmeas = false;
end

I = S.NextIx ~= 0;
Adj = sparse(S.SBIx(I),S.NextIx(I),1,nSBIx,nSBIx);


% nr of pixels in each subbasin
n = cell2mat(cellfun(@(x) numel(x),S.SBIdxList,'UniformOutput',false));
nagg = (speye(nSBIx)-Adj')\n;
if agg
    area = nagg*(cs^2);
else
    area = n*(cs^2);
end


if pixvalmeas
    A = A(:);
    
    % sum
    s = cell2mat(cellfun(@(x) sum(A(x(~isnan(A(x))))),S.SBIdxList,'UniformOutput',false));
    if agg
        s = (speye(nSBIx)-Adj')\s;
    end
    
    % mean
    m = cell2mat(cellfun(@(x) mean(A(x(~isnan(A(x))))),S.SBIdxList,'UniformOutput',false));
    if agg
        m = ((speye(nSBIx)-Adj')\(m.*n))./nagg;
    end
    
    % min
    mi = cell2mat(cellfun(@(x) min(A(x)),S.SBIdxList,'UniformOutput',false));
    if agg
        mi = cell2mat(cellfun(@(ix,ix2) min([mi(ix);mi(ix2)]),num2cell(S.SBIx),S.PreviousIx,'UniformOutput',false));
    end
    
    % max
    ma = cell2mat(cellfun(@(x) max(A(x)),S.SBIdxList,'UniformOutput',false));
    if agg
        ma = cell2mat(cellfun(@(ix,ix2) max([ma(ix);ma(ix2)]),num2cell(S.SBIx),S.PreviousIx,'UniformOutput',false));
    end
    
    r = ma-mi;
end

if ~pixvalmeas
    P = struct('Area',num2cell(area));
else
    P = struct('Area',num2cell(area),...
               'Mean',num2cell(m),...
               'Sum',num2cell(s),...
               'Min',num2cell(mi),...
               'Max',num2cell(ma),...
               'Range',num2cell(r)...          
           );
end
    







