function [DEM,bw] = ksdensity(DEM,x,y,varargin)

%KSDENSITY kernel density estimator for GRIDobj
%
% Syntax
%     
%     F = ksdensity(DEM,x,y)
%     F = ksdensity(DEM,x,y,pn,pv,...)
%
% Description
%
%     ksdensity returns a probability estimate f for the coordinates in x
%     and y, evaluated at the cell centers of the GRIDobj DEM.
%
%     See the function ksdensity (Statistics and Machine Learning Toolbox
%     for details).
%
% Input arguments
%
%     DEM     GRIDobj
%     x,y     coordinates
%     pn,pv   Parameter name, value pairs (see overloaded ksdensity function)
%             In addition, following values can be set
%             'useparallel'  {true} or false. If true, ksdensity will be
%                            run in parallel (requires Parallel Computing
%                            Toolbox)
%
% Output arguments
%
%     F       GRIDobj with probability estimates
%     bw      bandwidth
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     ix = randperm(prod(DEM.size),500);
%     [x,y] = ind2coord(DEM,ix);
%     [F,bw] = ksdensity(DEM,x,y,'bandwidth',1000);
%     imageschs(DEM,F)
%
% See also: GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 19. March, 2018


% go through varargin to find 'useparallel'
TF = strcmpi('useparallel',varargin);
ix = find(TF);
if ~isempty(ix)
    par = varargin{ix+1};
    varargin([ix ix+1]) = [];
else
    par = true;
end


[xc,yc] = getcoordinates(DEM);

if ~par
    [xc,yc] = meshgrid(xc,yc);
    [f,~,bw] = ksdensity([x,y],[xc(:) yc(:)],varargin{:});
    f = reshape(f,DEM.size);
else
    f = zeros(DEM.size);
    [xc,~] = meshgrid(xc,yc);
    bw = cell(size(f,2),1);
    parfor r = 1:size(f,2)
        [f(:,r),~,bw{r}] = ksdensity([x y],[xc(:,r) yc],varargin{:});
    end
    bw = bw{1};
end

DEM.Z = f;
DEM.name = 'kernel density';





