function L = laplacian(dem)

% construct a five-point laplacian matrix for solving diffusion equation


n        = numel(dem);
[ic,icd] = ixneighbors(dem,[],4);
L = sparse(ic,icd,1,n,n);
L = spdiags(sum(L,2),0,n,n) - L;

%% in order to calculate the final diffusion matrix, do following
% L = speye(n) + a * L;
% %where
% a = kappa*dt/(dx^2);



