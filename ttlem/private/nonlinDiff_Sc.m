function DEM= nonlinDiff_Sc(DEM,p,Sc,dt,forcingTerm,C)
% Calcualte nonlinear diffusion using the implicit method (Perron, 2011)
%
% Syntax
%
%       DEM_D= nonlinDiff_Sc(DEM,p,dt)
%
% Description
%
%       Calcualte nonlinear diffusion using the implicit method proposed by
%       Taylor Perron (2011): Numerical methods for nonlinear hillslope transport
%       laws. JGR-ES. 
%       No slopes can  exceed the 0.9 time the threshold slope to make sure
%       the algorithm is unconditionally stable. 
%       Source code can be found on: http://goo.gl/v1F0qs
%
% Input
%
%       DEM       DEM (digital elevation model) (GRIDobj)
%       p         parameter file
%       Sc          Threshold angle for landsliding in %
%       dt        time step 
%       forcingTerm         
%                 The forcing term is the combined impact of river incsion
%                 and uplift rates. 
%       C         River location
% 
% Output
%
%       DEM     DEM, updated for diffusion (GRIDobj)
%
% Example

% No slopes can exceed the threshold slope in the intital conditions 
[dx,dy]=gradient(DEM.Z,DEM.cellsize);
slope=sqrt(dx.^2+dy.^2);
slope(C.Z)=0;
if any(slope(:)>(Sc*0.9999))
    error('Slopes are exceeding 0.99 time Sc., choose another diffusion method!')
end

z=DEM.Z(:);

% boundaries
bdy=zeros(size(DEM.Z));
bdy(1,:)=1; % y boundaries
bdy(:,1)=1; % x boundaries
bdy(end,:)=1; % y boundaries
bdy(:,end)=1; % x boundaries
bdy(C.Z)=1;
bdy=bdy(:);

%% Following Perron:
N= prod(DEM.size);

% build LHS sparse matrix and RHS vector
[Ar, Ac, Av, i, B] = mexGetMats_BC_v1(z,bdy,dt,DEM.cellsize,DEM.cellsize,DEM.size(2),DEM.size(1),p.D,1/(Sc*Sc),forcingTerm);
A = sparse(Ar(1:i),Ac(1:i),Av(1:i),N,N);

% solve the system to get future elevations
% z = A\B;
[z,~] = pcg(A,B,p.DiffTol);
DEM.Z=reshape(z,DEM.size(1),DEM.size(2));
