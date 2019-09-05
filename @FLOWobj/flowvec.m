function [U,V] = flowvec(FD,L)

%FLOWVEC velocity vectors from FLOWobj
%
% Syntax
%
%     [U,V] = flowvec(FD)
%     [U,V] = flowvec(FD,L)
% 
% Description
%
%     flowvec returns the velocity components of the flow network in FD.
%     The components are normalized so that vector length is one. The
%     vector


U = GRIDobj(FD);
[X,Y] = getcoordinates(U);
X = single(X./FD.cellsize);
Y = single(Y./FD.cellsize);
[X,Y] = meshgrid(X,Y);

DX = X(FD.ixc) - X(FD.ix);
DY = Y(FD.ixc) - Y(FD.ix);

clear X Y

if ~ismulti(FD)
    %% Single flow directions
    V = U;
    U.Z(FD.ix) = DX;
    V.Z(FD.ix) = DY;
    
else
    %% Multiple flow directions    
    U.Z   = reshape(accumarray(FD.ix,DX,[prod(U.size) 1],...
                       @mean,nan('single')),U.size);
    V     = U;
    V.Z   = reshape(accumarray(FD.ix,DY,[prod(U.size) 1],...
                       @mean,nan('single')),U.size);  
end

% Normalize vectors to have unit length
LL = sqrt(U.^2 + V.^2);
U  = U./LL;
V  = V./LL;

if nargin == 2
    U = U.*L;
    V = V.*L;
end
 