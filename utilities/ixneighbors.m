function [ic,icd] = ixneighbors(varargin)

% neighbor indexing of a matrix
%
% Syntax
%
%     [ic,icd] = ixneighbors(A)
%     [ic,icd] = ixneighbors(A,ix)
%     [ic,icd] = ixneighbors(A,I)
%     [ic,icd] = ixneighbors(...,conn)
%
% Description
%
%     ixneighbors returns the indices of neighbor cells in a n*m matrix A.
%     ic and icd are column vectors of same length where ic contains the
%     indices of cells in A and icd contains the indices of the neighbor
%     cells.
%
%     ixneighbors(A)      returns all neighbors of all cells in A
%     ixneighbors(A,ix)   returns all neighbors of the cells in index 
%                         vector ix
%     ixneighbors(A,I)    returns all neighbors of the cells in the logical
%                         matrix I that are TRUE. I must be same size as A
%     ixneighbors(...,conn) lets you specify the type of connectivity. conn
%                         is either 4 or 8. By default, ixneighbors uses
%                         an eight-connectivity
%     
%
%     ixneighbors handles NaNs. Hence, it discards cells in A that are NaN 
%     both in ic and icd.
%
% Example 1
% 
%     A = magic(4);
%     A(2,2) = NaN
% 
%     A =
% 
%         16     2     3    13
%          5   NaN    10     8
%          9     7     6    12
%          4    14    15     1
% 
%     [ic,icd] = ixneighbors(A,3)
% 
%     ic =
% 
%          3
%          3
%          3
%          3
% 
% 
%     icd =
% 
%          7
%          4
%          2
%          8
%
%
% Example 2
%
%     % Construct a sparse adjacency matrix S
% 
%     A = peaks(100);
%     A(A<0) = NaN;
%     nrc = numel(A);
%     [ic,icd] = ixneighbors(A);
%     S = sparse(ic,icd,ones(numel(icd),1),nrc,nrc);
%     spy(S)
%
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009



% handle input and error checking
if nargout~=2;
    error('wrong number of output arguments')
end

X   = varargin{1};
siz = size(X);
nrc = numel(X);

In  = isnan(X);


if nargin==1;   
    method = 'getall';
    nhood  = 8;
elseif nargin==2 || nargin==3;
    ix = varargin{2};
    if isempty(ix);
        method = 'getall';
    else
        method = 'getsome';
        if islogical(ix)
            if size(X) ~= size(ix);
                error('if I is logical I and X must have same size')
            end
        else
            ixvec = ix(:);
            ix = false(siz);
            ix(ixvec) = true;
        end
        ix = ~In & ix;
    end
    
    
    if nargin==3;
        nhood = varargin{3};
        if ~ismember(nhood(1),[4 8]);
            error('conn must be either 4 or 8')
        else
            nhood = nhood(1);
        end
    else
        nhood = 8;
    end
    
else 
    error('wrong number of input arguments')
end
 
% replace values in X by index vector

X     = reshape((1:nrc)',siz);
X(In) = NaN;
% Pad array
ic  = nan(siz(1)+2,siz(2)+2);
ic(2:end-1,2:end-1) = X;


switch method
    case 'getall'
        I   = ~isnan(ic);
    case 'getsome'
        % Pad logical array
        I = false(siz(1)+2,siz(2)+2);
        I(2:end-1,2:end-1) = ix;
end

icd = zeros(nnz(I),nhood);


% Shift logical matrix I across the neighbors
icd(:,1) = ic(I(:,[end 1:end-1]));                % shift to the right                    
icd(:,2) = ic(I([end 1:end-1],:));                % shift down       
icd(:,3) = ic(I(:,[2:end 1]));                    % shift left        
icd(:,4) = ic(I([2:end 1],:));                    % shift up      

if nhood==8;

icd(:,5) = ic(I([2:end 1],[end 1:end-1]));        % shift up and right        
icd(:,6) = ic(I([2:end 1],[2:end 1]));            % shift up and left        
icd(:,7) = ic(I([end 1:end-1],[end 1:end-1]));    % shift down and right        
icd(:,8) = ic(I([end 1:end-1],[2:end 1]));        % shift down and left

end
% Create output
ic = repmat(ic(I(:)),nhood,1);
icd = icd(:);

% Remove NaNs in neighbors
i = isnan(icd);

ic(i) = [];
icd(i) = [];
               
        