function M = FLOWobj2M(FD)

% convert instance of FLOWobj to flow direction matrix 
%
% Syntax
%
%     M = FLOWobj2M(FD);
%
% Description
%
%     FLOWobj2M converts an instance of FLOWobj to the flow direction
%     matrix M as used in TopoToolbox 1.
%
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. February, 2013


nrc = prod(FD.size);
switch FD.type
    case 'multi'
        M = sparse(double(FD.ix),double(FD.ixc),FD.fraction,nrc,nrc);
    case 'single'
        M = sparse(double(FD.ix),double(FD.ixc),1,nrc,nrc);
end
end