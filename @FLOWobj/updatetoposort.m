function FD = updatetoposort(FD)

%UPDATETOPOSORT update topological sorting
%
% Syntax
%
%     FDu = updatetoposort(FD)
%
% Description
%
%     updatetoposort topologically sorts edges in a FLOWobj. This may be 
%     necessary if edge vectors were manipulated.
%
% Input arguments
%
%     FD   FLOWobj
%
% Output arguments
%
%     FDu  updated FLOWobj
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 21. March, 2017

switch FD.type
    case 'single'
        nr= prod(FD.size);
        D = digraph(FD.ix,FD.ixc);
        p = toposort(D);
        
        IX = zeros(nr,1,'uint32');
        IX(FD.ix) = FD.ix;
        IXC = zeros(nr,1,'uint32');
        IXC(FD.ix) = FD.ixc;

        p   = p(:);
        IX  = IX(p);
        IXC = IXC(p);
        IX  = nonzeros(IX);
        IXC = nonzeros(IXC);
        FD.ix  = uint32(IX(:));
        FD.ixc = uint32(IXC(:));
        
    case 'multi'
        
        M = FLOWobj2M(FD);
        G = digraph(M);
        p = toposort(G);
        [FD.ix,FD.ixc,FD.fraction] = find(M(p,p));
        
        p = p(:);
        FD.ixc = uint32(p(FD.ixc));
        FD.ix  = uint32(p(FD.ix));
        
end



