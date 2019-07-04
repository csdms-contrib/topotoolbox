function ind=getIndices(A, BC)


ind.inner=false(size(A));
ind.inner((1+BC.nbGhost):(end-BC.nbGhost),(1+BC.nbGhost):(end-BC.nbGhost)) = true;

if BC.nbGhost==1
    ind.left=false(size(A));
    ind.left(2:end-1,1:end-2) = true;
    ind.right=false(size(A));
    ind.right(2:end-1,3:end) = true;
    ind.top=false(size(A));
    ind.top(1:end-2,2:end-1) = true;
    ind.bottom=false(size(A));
    ind.bottom(3:end,2:end-1) = true;
elseif BC.nbGhost==2
    ind.left=false(size(A));
    ind.left(3:end-2,2:end-3) = true;
    ind.right=false(size(A));
    ind.right(3:end-2,4:end-1) = true;
    ind.top=false(size(A));
    ind.top(2:end-3,3:end-2) = true;
    ind.bottom=false(size(A));
    ind.bottom(4:end-1,3:end-2) = true;
    
    ind.left2=false(size(A));
    ind.left2(3:end-2,1:end-4) = true;
    ind.right2=false(size(A));
    ind.right2(3:end-2,5:end) = true;
    ind.top2=false(size(A));
    ind.top2(1:end-4,3:end-2) = true;
    ind.bottom2=false(size(A));
    ind.bottom2(5:end,3:end-2) = true;
    
end
end