function D = flowdistanceds(M,X,Y,dem,ix)

% compute maximum downstream flow distance from each cell in the DEM
% 
% Syntax
%
%     D = flowdistanceds(M,X,Y)
%     D = flowdistanceds(M,X,Y,dem)
%     D = flowdistanceds(M,X,Y,ix)
%     D = flowdistanceds(M,X,Y,dem,ix)
%
% Description
%
%     flowdistance computes the maximum downstream flow path length from 
%     each cell to the catchment outlet. The computed distance is the
%     euclidean 2D or 3D distance along the flowpath indicated in the
%     single flow direction matrix M. Note that for correct distances
%     the units of the coordinate matrices and the elevation model must be 
%     the same.
%
% Input
%
%     M       single flow direction matrix
%     X, Y    coordinate matrices as generated by meshgrid
%     dem     digital elevation model (when supplied the euclidean 3D 
%             distance is calculated)
%     ix      linear index of the cell from which the distance shall be
%             computed. ix must be a nonnegative scalar integer. If ix is
%             not supplied, flowdistance calculates the distance to each
%             catchment outlet.
%
% Output
%
%     D       distance matrix
%
% Example
%
%     load exampleDEM
%     M = flowdir_single(dem);
%     ix = 10365; % location of interest
%     D = flowdistanceds(M,X,Y,[],ix);
%     imagesc(X(1,:),Y(:,1),D);axis image; axis xy
%     colorbar
%     title('distance [m] from yellow circle downstream')
%     hold on
%     plot(X(ix),Y(ix),'yo')
%     hold off
%     
%
%
% See also: FLOWDIR, DEPENDENCEMAP
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009



siz = size(X);
nrc = numel(X);

% check input arguments
if nargin <= 2;
    error('TopoToolbox:incorrectinput',...
              'Too few input arguments')
elseif nargin == 3;
    algo = 'all';
    flag3d = false;
elseif nargin == 4;
    if isscalar(dem);
        algo = 'one';
        ix   = dem;
        flag3d = false;
    else
        algo = 'all';
        flag3d = true;    
    end
    
elseif nargin == 5;
    algo = 'one';
    if isempty(dem);
        flag3d = false;
    else
        flag3d = true;
    end
    if ~isscalar(ix);
        error('TopoToolbox:incorrectinput',...
              'ix must be a scalar')
    end
else
    error('TopoToolbox:incorrectinput',...
          'wrong number of input arguments')
end

% do the matrices have same size
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

% does the flow direction matrix correspond to the DEM
if nrc~=size(M,1) || nrc~=size(M,2);
    error('TopoToolbox:incorrectinput',...
          'M must be a square, sparse matrix with size [numel(X) numel(X)]')
end

% check if single flow direction Matrix
% is used
if any(sum(spones(M),2)>1);
    error('TopoToolbox:incorrectinput',...
          'single flow direction matrix must be used')
end


switch algo
    case 'one'
        
        % get linear indices of flow connectivity
        [ic,icd] = find(M);
        % calculate distance between nodes
        f = hypot(X(ic)-X(icd),Y(ic)-Y(icd));
        
        % 3d?
        if flag3d
            f = hypot(f,dem(ic)-dem(icd));
        end
        
        
        B = zeros(nrc,1);
        B(ic) = f;
        
        % set distances to zero that are not upslope to ix
        
        C = zeros(nrc,1);
        
        if ix == 0 || ix >nrc;
            error('TopoToolbox:incorrectinput',...
                ['ix must range between 1 and ' num2str(nrc)])
        end
        
        C(ix) = 1;
        I = (speye(nrc)-M')\C;
        B(~(I>0)) = 0;
        B(ix) = 0;
        
        D = (speye(nrc)-M')\B;
        D = reshape(D,siz);
        
        
        
        
        
    case 'all'
               
        % defaults
        cs  = abs(Y(1)-Y(2));
        
        % calculate upstream flowdistance
        if flag3d;
            D = flowdistance(M,X,Y,dem);
        else
            D = flowdistance(M,X,Y);
        end
        
        % Calculation of Slope Length
        [ic,icd]   = find(M);
        [dummy,ix] = sort(D(ic),'descend');
        
        clear dummy
        
        ic          = ic(ix);
        icd         = icd(ix);
        [icdd,icdd] = ismember(icd,ic);
        
        
        % initialize length-slope
        SL  = zeros(siz);
        % Find non-receivers (cells on ridges)
        NG  = sum(M,1)' == 0;
        % slope length in non-receiver cells in half the
        % cellsize
        SL(NG) = cs/2;
        
        NG  = NG(ic);
        NG  = find(NG);
        % nr of nonreceivers
        nNG = numel(NG);
        
        NGix = 1;
        
        D = abs(D(ic)-D(icd));
        
        IX  = 1;
        
        
        % not vectorized
        while NGix <= nNG;
            
            % distance
            ln = D(IX) + SL(ic(IX));
            
            if ln > SL(icd(IX));
                SL(icd(IX)) = ln;
                
                % next row to work on
                IX = icdd(IX);
                
                flaggoon = IX~= 0;
            else
                flaggoon = false;
            end
            
            if ~flaggoon;
                
                NGix = NGix+1;
                
                if NGix>nNG;
                else
                    IX  = NG(NGix);
                end
            end
        end
        
        D = SL;
end


