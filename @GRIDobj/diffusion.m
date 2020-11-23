function DEM = diffusion(DEM,varargin)

%DIFFUSION Solve the diffusion equation
%
% Syntax
%
%     DEMd = diffusion(DEM,D,dt)
%
% Description
%
%     This function uses an implicit scheme to solve the linear diffusion
%     equation for a DEM.
%
% Input arguments
%
%     DEM     digital elevation model
%
%     Parameter name/value pairs
%
%     D           diffusivity (m^2 /y) (default = 1)
%     timespan    duration (y) (default = 1000)
%     numsteps    number of iterations (default = 5)
%     streamnet   STREAMobj
%     solver      'pcg' (default) or '\'
%     pcgtol      1e-6 (default)
%     uplift      GRIDobj or scalar (default = 0)
%
%
% Output arguments
%
%     DEMd    diffused DEM
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMd = DEM;
%     for r = 1:10;
%     DEMd = diffusion(DEMd); 
%     imageschs(DEMd); 
%     drawnow; 
%     end
%     figure
%     imageschs(DEM,DEMd-DEM)
%
% See also: GRIDobj/filter, ttlem
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 20. October, 2020

p = inputParser;
addParameter(p,'timespan',1000)
addParameter(p,'numsteps',5)
addParameter(p,'D',1)
addParameter(p,'streamnet',[])
addParameter(p,'solver','pcg')
addParameter(p,'pcgtol',1e-6)
addParameter(p,'uplift',0)
parse(p,varargin{:});

D = p.Results.D;
numsteps = p.Results.numsteps;
dt = p.Results.timespan / numsteps;



nrc = prod(DEM.size);
% get neighbor indices
[ic,icd] = ixneighbors(DEM.Z,[],4);

% remove indices to nan-cells
I = isnan(DEM.Z(ic)) | isnan(DEM.Z(icd));
ic(I) = [];
icd(I) = [];

% calculate laplacian
L        = sparse(ic,icd,1,nrc,nrc);
L        = spdiags(sum(L,2),0,nrc,nrc) - L;
% and the diffusion matrix
if isa(D,'GRIDobj')
	D = D.Z(:);
else
	D = repmat(D,numel(DEM.Z),1);
end

if isempty(p.Results.streamnet)
	D  = speye(nrc) + spdiags(D,0,nrc,nrc)*dt/(2*DEM.cellsize^2)*L;
else
	S  = +STREAMobj2GRIDobj(p.Results.streamnet);
	D  = speye(nrc) + spdiags(D,0,nrc,nrc)*dt/(2*DEM.cellsize^2)*spdiags(1-S.Z(:),0,nrc,nrc)*L;
end
% must work with doubles
% remember class
c        = class(DEM.Z);
% ensure double
DEM.Z    = double(DEM.Z);

% set nan values to zero
I        = isnan(DEM.Z);
DEM.Z(I) = 0;

% solve
Z1 = DEM.Z(:);

if isa(p.Results.uplift,'GRIDobj')
    u = p.Results.uplift.Z;
else 
    u = p.Results.uplift;
    u = repmat(u,DEM.size);
end

if ~isempty(p.Results.streamnet)
    u(S.Z>0) = 0;
end
u = u(:);
u = u/1000 * dt;

for r = 1:numsteps
    switch p.Results.solver
        case '\'
            Z1 = D\(Z1+u);
        case 'pcg'
            [Z1,~] = pcg(D,Z1+u,p.Results.pcgtol,[],[],[],Z1);
        otherwise
            error('unknown solver')
    end
end
% reshape
DEM.Z = reshape(Z1,DEM.size);
% reset values to nan
DEM.Z(I) = nan;
% reset input class
DEM.Z = cast(DEM.Z,c);