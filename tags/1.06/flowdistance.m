function [D,L] = flowdistance(M,X,Y,dem,ix)

% compute flow distance to cell or to catchment outlet
% 
% Syntax
%
%     [D,L] = flowdistance(M,X,Y)
%     [D,L] = flowdistance(M,X,Y,dem)
%     [D,L] = flowdistance(...,ix)
%     [D,L] = flowdistance(...,I)
%
% Description
%
%     flowdistance computes the flow path length from all cells to the
%     catchment outlet or specified cells. The computed distance is the
%     euclidean 2D or 3D distance along the flowpath represented by the 
%     single flow direction matrix M. Note that for correct distances 
%     the units of the coordinate matrices and the elevation model must
%     be the same. The second output of flowdistance contains the linear
%     indices of the closest catchment outlet (in case of three or four
%     input arguments) or the closest nonzero values in I (or ix).
%
%     D can be regarded as a special case of a distance transform of the
%     DEM based a quasi-euclidean metric along flowpaths (see bwdist for 
%     further information on distance transforms.
%
%     Assuming constant flow velocity in the catchment isochores can be
%     calculated from D.
%
% Input
%
%     M       single flow direction matrix
%     X, Y    coordinate matrices as generated by meshgrid
%     dem     digital elevation model (when supplied the euclidean 3D 
%             distance is calculated). Supply an empty matrix ([]) if you 
%             wish to calculate 2D distances.
%     ix      linear index of the cell from which the distance shall be
%             computed. ix must be a nonnegative scalar integer. 
%     I       logical matrix same size as dem. flowdistance than computes
%             distance along the flowpath from zeros cells in I to the
%             nearest nonzero cell.
%
%             If ix or I are  not supplied, flowdistance calculates
%             the distance to each catchment outlet.
%
% Output
%
%     D       distance matrix
%     L       L contains the linear index of the nearest nonzero pixel of
%             I.
%
% Example 1
%
%     load exampleDEM
%     M = flowdir_single(dem);
%     D = flowdistance(M,X,Y,dem);
%     imagesc(X(1,:),Y(:,1),D);axis image; axis xy
%     title('flow path distance [m] to grid edges')
%     colorbar
%
%
% Example 2
%
%     % create a simplified hydrograph, where an
%     % amount of one is in each cell and is routed downstream
%     % with steady velocity
%
%     load exampleDEM
%     M = flowdir_single(fillsinks(dem));
%     A = flowacc(M,size(dem));
%     [ix,ix] = max(A(:));
%     D = flowdistance(M,X,Y,dem,ix);
%     imagesc(X(1,:),Y(:,1),D);axis image; axis xy
%     figure
%     hist(D(D~=0),500);
%     xlabel('time*velocity')
%     ylabel('discharge')
%
% Example 3
%
%     % investigate the relation between gradients of hillslopes 
%     % adjacent to the drainage network (within a distance of 300 m
%     % along the flowpath)
% 
%     load exampleDEM
%     M = flowdir_single(fillsinks(dem));
%     A = flowacc(M,size(dem));
%     G = gradient8(dem,abs(Y(1)-Y(2)));
%     D = flowdistance(M,X,Y,dem,A>100);
%     I = D>0 & D <=sqrt((2*(Y(1)-Y(2)))^2);
%     semilogx(A(I),G(I),'.')
%
% See also: FLOWDIR, DEPENDENCEMAP, BWDIST, FLOWPATHBUFFER
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 5. February, 2010


error(nargchk(2, 5, nargin))


siz = size(X);
nrc = numel(X);

% check input arguments
if nargin == 2;
    siz = X;
    nrc = prod(siz);
    flag3d = false;
    flagix = false;
    
elseif nargin == 3;
    flagix = false;
    flag3d = false;
else
    flag3d = ~isempty(dem);
    flagix = false;
    
    if nargin == 5
        flagix = true;
    end
end
    
    
% do the matrices have same size
if nargin > 2;
    
    if ~flag3d
        if ~isequal(siz,size(Y))
            error('TopoToolbox:incorrectinput',...
                'X and Y must have same size')
        end
    else
        if ~isequal(siz,size(Y),size(dem))
            error('TopoToolbox:incorrectinput',...
                'X, Y and dem must have same size')
        end
    end
end
% does the flow direction matrix correspond to the DEM
if nrc~=size(M,1) || nrc~=size(M,2);
    error('TopoToolbox:incorrectinput',...
          'M must be a square, sparse matrix with size [numel(X) numel(X)]')
end

% check if single flow direction Matrix
% is used
M = spones(M);
if any(sum(M,2)>1);
    error('TopoToolbox:incorrectinput',...
          'single flow direction matrix must be used')
end


% check if a logical matrix same size as X is supplied as fifth argument
if flagix
    if isequal(size(ix),siz)
        II = ix>0;
    else
        II = false(siz);
        II(ix) = true;
    end
    
    % disconnect values where B == true from downward neighbors
    M = spdiags(~II(:),0,nrc,nrc)*M;
    % calculate dependence map
    I = dependencemap(M,II);
    % remove connectivity in downslope areas that are not included in the
    % calculation
    M = spdiags(I(:),0,nrc,nrc)*M;
    
end
M = spones(M);

if nargin > 2;
    % get linear indices of flow connectivity
    [ic,icd] = find(M);
    % calculate distance between nodes
    f = hypot(X(ic)-X(icd),Y(ic)-Y(icd));
    if flag3d
        f = hypot(f,dem(ic)-dem(icd));
    end
    
    %
    B     = zeros(nrc,1);
    B(ic) = f;
else
    B = ones(nrc,1);
end

D = (speye(nrc)-M)\B;
D = reshape(D,siz);

if nargout == 2;
    % identify cells that only receive
    I = (sum(M,1) > 0)' & sum(M,2) == 0;
    L = zeros(nrc,1);
    L(I) = find(I);
    
    L = (speye(nrc)-M)\L;
    L = reshape(L,siz);
    L(II) = false;
end



