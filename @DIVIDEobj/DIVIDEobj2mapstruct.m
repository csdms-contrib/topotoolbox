function MS = DIVIDEobj2mapstruct(D,DEM,seglen,varargin)
%DIVIDEPROPS obtain divide properties from GRIDobj
%
% Syntax
%
%     D = DIVIDEobj2mapstruct(D,DEM,seglength)
%     D = DIVIDEobj2mapstruct(D,DEM,seglength,...
%             {'fieldname1' var1 aggfunction1},...
%             {'fieldname2' var2 aggfunction2})
%
% Description
%
%     DIVIDEobj2mapstruct creates a geographic data structure that can be
%     exported as a shapefile. The divide segments are given attributes
%     that are derived from the DEM, and optionally other GRIDobjs or 
%     lists of values associated to the divide edges.
%     A segment length value of zero (0) yields for each feature the
%     original divide segments. For all other (positive) segment length 
%     values, the divide segments will be chopped into 
%
%
% Input
%
%     D         instance of DIVIDEobj
%     DEM       digital elevation model (GRIDobj)
%     seglength length of divide segments
%
%     {'fieldname' var aggfunction}
%               the triplets fieldname (char), var (GRIDobj), and
%               aggfunction ('min','max','mean','diff') 
%     
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     FA = flowacc(FD);
%     ST = STREAMobj(FD,FA>500);
%     D = DIVIDEobj(FD,ST);
%     D = divorder(D,'topo');
%     DX = flowdistance(FD,ST);
%     DZ = vertdistance2stream(FD,ST,DEM);
%     DZ.Z(DZ.Z<0) = 0;
%     DZ.Z(isinf(DZ.Z)) = nan;
%     % create mapping structure
%     MS = DIVIDEobj2mapstruct(D,DEM,1000,...
%         {'hr_mean' DZ 'mean'},{'hr_diff' DZ 'diff'},...
%         {'fdist_mean' DX 'mean'},{'fdist_diff' DX 'diff'});
%     for i = 1 : length(MS)
%         MS(i).dai = MS(i).hr_diff./MS(i).hr_mean./2;
%     end
%     % visualize
%     do = [MS.order];
%     dai = [MS.dai];
%     symbolspec = makesymbolspec('Line',...
%         {'order',[2 max(do)],'Linewidth',[0.5 6]},...
%         {'dai',[0 1],'Color',flip(hot)});
%     figure
%     imageschs(DEM,DEM)
%     hc = colorbar;
%     hc.Label.String = 'Elevation (m)';
%     hold on
%     ix = do>1 & not(isnan(dai));
%     mapshow(MS(ix),'SymbolSpec',symbolspec);
%     title('Divide asymmetry index: white=low -> red=high')
%     % export as shapefile, if needed: 
%     % shapewrite(MS,'big_tujunga_drainage_divides.shp')
%    
%
% See also: FLOWobj, divides
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: April 2020


if not(D.issorted)
    error('Divide object not sorted. Use function SORT first.')
end

if isempty(D.distance)
    warning('Divide object has no distance. Comnputing distance now.')
    D = divdist(D);
end
if isempty(D.distance)
    warning('Divide object has no ordertype. Comnputing Topo-order now.')
    D = divorder(D,'topo');
end


% Get vectors 
[x,y] = ind2coord(D,vertcat(D.IX));
do = D.order;
dist = D.distance;

% Adjust segment length 
if seglen>0 && not(isnan(seglen))
    % (copied from function STREAMobj2mapstruct)
    sep = isnan(x);
    d  = zeros(size(x));
    for r = 2:numel(d) 
        if sep(r-1)
            d(r) = 0;
        elseif sep(r)
            d(r) = nan;
        else
            d(r) = d(r-1)+sqrt((x(r)-x(r-1)).^2 + (y(r)-y(r-1)).^2);
        end
    end
    % try to get equal length segments being equally distributed over an
    % existing segment.
    % Use equal quantiles
    ixsep = find(sep);
    nrsep = numel(ixsep);
    ixsep = [0;ixsep];
    newsep   = zeros(size(d));
    for r = 1:nrsep
        dd = d(ixsep(r)+1 : ixsep(r+1)-1);
        ddmax = max(dd);
        nrnewsegs = max(round(ddmax/seglen),1);
        if nrnewsegs > 1
            [~,segtemp] = histc(dd,...
                linspace(0,ddmax + DEM.cellsize*.01,nrnewsegs));
            segtemp = [0;diff(segtemp)];
            newsep(ixsep(r)+1 : ixsep(r+1)-1) = segtemp;
        end
    end
    % now set new separators according to newsep
    ind = find(newsep > 0);
    x  = insertrows(x,x(ind),ind);
    y  = insertrows(y,y(ind),ind);
    do = insertrows(do,do(ind),ind);
    dist = insertrows(dist,dist(ind),ind);
    newsep = insertrows(newsep,0,ind);
    % insert nans
    ind = find(newsep > 0);
    x  = insertrows(x,nan,ind);
    y  = insertrows(y,nan,ind);
    do = insertrows(do,nan,ind);
    dist = insertrows(dist,nan,ind);
end


% Calculate across divide attributes 
nvar = length(varargin);
VAR = cell(1,3+nvar);
varname = cell(1,3+nvar);

% Find pixels on either side of ridgeline
x1 = [NaN; x(1:end-1)];
x2 = [NaN; x(2:end)];
y1 = [NaN; y(1:end-1)];
y2 = [NaN; y(2:end)];
dx = x1-x2;
dy = y1-y2;
hcs = DEM.cellsize/2;
ix = dx==0; % vertical link
iy = dy==0; % horizontal link
meanx = (x1+x2)./2;
meany = (y1+y2)./2;
px = meanx + hcs.*ix;
qx = meanx - hcs.*ix;
py = meany + hcs.*iy;
qy = meany - hcs.*iy;
pix = coord2ind(DEM,px,py);
qix = coord2ind(DEM,qx,qy);

% allocate space
px1 = nan(size(pix));
px2 = px1;
nx = ~isnan(pix) & ~isnan(qix);

% get elevation
px1(nx) = DEM.Z(pix(nx));
px2(nx) = DEM.Z(qix(nx));
minpx = min([px1,px2],[],2);
maxpx = max([px1,px2],[],2);
VAR{1} = min([minpx maxpx],[],2);
varname{1} = 'z_min';
VAR{2} = max([minpx maxpx],[],2);
varname{2} = 'z_max';
VAR{3} = mean([minpx maxpx],2);
varname{3} = 'z_mean';
ct = 3;


% get other grid values
if ~isempty(varargin)
    for i = 1 : nvar
        GRID = varargin{i}{2};
        if isa(GRID,'GRIDobj')
            
            px1(nx) = GRID.Z(pix(nx));
            px2(nx) = GRID.Z(qix(nx));
            switch varargin{i}{3}
                case 'mean'
                    v = nanmean([px1,px2],2);
                case 'max'
                    v = nanmax([px1,px2],[],2);
                case 'min'
                    v = nanmin([px1,px2],[],2);
                case 'diff'
                    v = abs(diff([px1,px2],1,2));
            end
        else
            if length(GRID)==length(D.IX)
                v = GRID;
            else
                error('Attribute variable size does not match DIVIDEobj length');
            end
        end
        ct = ct+1;
        VAR{ct} = v;
        varname{ct} = varargin{i}{1};
    end
end
M = cell2mat(VAR);


% Create the mapping structure 
ix      = find(isnan(x));
nrlines = numel(ix);

MS = struct('Geometry','Line',...
    'X',cell(nrlines,1),...
    'Y',cell(nrlines,1));

IXs = [1;ix(1:end-1)+1];
IXe = ix-1;

for r = 1:nrlines
    MS(r).ID = r;
    MS(r).X = [x(IXs(r):IXe(r));NaN];
    MS(r).Y = [y(IXs(r):IXe(r));NaN];
    if not(isempty(do))
        MS(r).order = unique(do(IXs(r):IXe(r))); % if not unique, something wrong
    end
    if not(isempty(dist))
        MS(r).distance = mean(dist(IXs(r):IXe(r)));
    end
    MS(r).length = max(getdistance(MS(r).X,MS(r).Y));
    % calculate orientation
    dx = MS(r).X(end-1)-MS(r).X(1);
    dy = MS(r).Y(end-1)-MS(r).Y(1);
    alpha = atand(dx./dy);
    if alpha<0; alpha = alpha+180; end
    MS(r).azimuth = alpha;
    for k = 1 : ct
        MS(r).(varname{k}) = double(nanmean(M(IXs(r):IXe(r),k)));
    end
end

end % main function


%% INSERTROWS by JOS 
% www.mathworks.com/matlabcentral/fileexchange/9984-insertrows-v2-0-may-2008
function [C,RA,RB] = insertrows(A,B,ind)
% INSERTROWS - Insert rows into a matrix at specific locations
%   C = INSERTROWS(A,B,IND) inserts the rows of matrix B into the matrix A at
%   the positions IND. Row k of matrix B will be inserted after position IND(k)
%   in the matrix A. If A is a N-by-X matrix and B is a M-by-X matrix, C will
%   be a (N+M)-by-X matrix. IND can contain non-integers.
%
%   If B is a 1-by-N matrix, B will be inserted for each insertion position
%   specified by IND. If IND is a single value, the whole matrix B will be
%   inserted at that position. If B is a single value, B is expanded to a row
%   vector. In all other cases, the number of elements in IND should be equal to
%   the number of rows in B, and the number of columns, planes etc should be the
%   same for both matrices A and B. 
%
%   Values of IND smaller than one will cause the corresponding rows to be
%   inserted in front of A. C = INSERTROWS(A,B) will simply append B to A.
%
%   If any of the inputs are empty, C will return A. If A is sparse, C will
%   be sparse as well. 
%
%   [C, RA, RB] = INSERTROWS(...) will return the row indices RA and RB for
%   which C corresponds to the rows of either A and B.
%
%   Examples:
%     % the size of A,B, and IND all match
%        C = insertrows(rand(5,2),zeros(2,2),[1.5 3]) 
%     % the row vector B is inserted twice
%        C = insertrows(ones(4,3),1:3,[1 Inf]) 
%     % matrix B is expanded to a row vector and inserted twice (as in 2)
%        C = insertrows(ones(5,3),999,[2 4])
%     % the whole matrix B is inserted once
%        C = insertrows(ones(5,3),zeros(2,3),2)
%     % additional output arguments
%        [c,ra,rb] = insertrows([1:4].',99,[0 3]) 
%        c.'     % -> [99 1 2 3 99 4] 
%        c(ra).' % -> [1 2 3 4] 
%        c(rb).' % -> [99 99] 
%
%   Using permute (or transpose) INSERTROWS can easily function to insert
%   columns, planes, etc:
%
%     % inserting columns, by using the transpose operator:
%        A = zeros(2,3) ; B = ones(2,4) ;
%        c = insertrows(A.', B.',[0 2 3 3]).'  % insert columns
%     % inserting other dimensions, by using permute:
%        A = ones(4,3,3) ; B = zeros(4,3,1) ; 
%        % set the dimension on which to operate in front
%        C = insertrows(permute(A,[3 1 2]), permute(B,[3 1 2]),1) ;
%        C = ipermute(C,[3 1 2]) 
%
%  See also HORZCAT, RESHAPE, CAT

% for Matlab R13
% version 2.0 (may 2008)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History:
% 1.0, feb 2006 - created
% 2.0, may 2008 - incorporated some improvements after being selected as
% "Pick of the Week" by Jiro Doke, and reviews by Tim Davis & Brett:
%  - horizontal concatenation when two arguments are provided
%  - added example of how to insert columns
%  - mention behavior of sparse inputs
%  - changed "if nargout" to "if nargout>1" so that additional outputs are
%    only calculated when requested for

narginchk(2,3);

if nargin==2,
    % just horizontal concatenation, suggested by Tim Davis
    ind = size(A,1) ;
end

% shortcut when any of the inputs are empty
if isempty(B) || isempty(ind),    
    C = A ;     
    if nargout > 1,
        RA = 1:size(A,1) ;
        RB = [] ;
    end
    return
end

sa = size(A) ;

% match the sizes of A, B
if numel(B)==1,
    % B has a single argument, expand to match A
    sb = [1 sa(2:end)] ;
    B = repmat(B,sb) ;
else
    % otherwise check for dimension errors
    if ndims(A) ~= ndims(B),
        error('insertrows:DimensionMismatch', ...
            'Both input matrices should have the same number of dimensions.') ;
    end
    sb = size(B) ;
    if ~all(sa(2:end) == sb(2:end)),
        error('insertrows:DimensionMismatch', ...
            'Both input matrices should have the same number of columns (and planes, etc).') ;
    end
end

ind = ind(:) ; % make as row vector
ni = length(ind) ;

% match the sizes of B and IND
if ni ~= sb(1),
    if ni==1 && sb(1) > 1,
        % expand IND
        ind = repmat(ind,sb(1),1) ;
    elseif (ni > 1) && (sb(1)==1),
        % expand B
        B = repmat(B,ni,1) ;
    else
        error('insertrows:InputMismatch',...
            'The number of rows to insert should equal the number of insertion positions.') ;
    end
end

sb = size(B) ;

% the actual work
% 1. concatenate matrices
C = [A ; B] ;
% 2. sort the respective indices, the first output of sort is ignored (by
% giving it the same name as the second output, one avoids an extra 
% large variable in memory)
[~,abi] = sort([[1:sa(1)].' ; ind(:)]) ;
% 3. reshuffle the large matrix
C = C(abi,:) ;
% 4. reshape as A for nd matrices (nd>2)
if ismatrix(A),
    sc = sa ;
    sc(1) = sc(1)+sb(1) ;
    C = reshape(C,sc) ;
end

if nargout > 1,
    % additional outputs required
    R = [zeros(sa(1),1) ; ones(sb(1),1)] ;
    R = R(abi) ;
    RA = find(R==0) ;
    RB = find(R==1) ;
end
end



