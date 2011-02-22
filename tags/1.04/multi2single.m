function M = multi2single(M)

% convert multiple to single flow direction matrix
%
% Syntax
% 
%     Ms = multi2single(M)
%
% Description
%
%     Various functions in the TopoToolbox only work with single flow
%     direction matrices. multi2single converts a multiple flow direction
%     matrix to a single flow direction matrix.
%
% Input
%
%     M         multiple flow direction matrix
%
% Output
%
%     Ms        single flow direction matrix
%
% Example
%
%     load exampleDEM
%     [A,M] = ezflowacc(X,Y,dem);
%     subplot(1,2,1)
%     imagesc(X(1,:),Y(:,2),A); axis image; axis xy
%     title('multiple flow direction')
% 
%     Ms = multi2single(M);
%     As = flowacc(Ms,size(dem));
%     subplot(1,2,2)
%     imagesc(X(1,:),Y(:,2),As); axis image; axis xy
%     title('single flow direction')
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009


% check input
error(nargchk(1, 1, nargin));

nrc     = size(M,1);
% identify largest element in each row of M
[m,icd] = max(M,[],2);

I       = m~=0;

ic      = (1:nrc)';

M = sparse(ic(I),icd(I),1,nrc,nrc);



