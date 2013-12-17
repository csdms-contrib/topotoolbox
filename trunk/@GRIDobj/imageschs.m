function rgb = imageschs(DEM,A,varargin)

% plot hillshade image with overlay
%
% Syntax
%
%     imageschs(DEM)
%     imageschs(DEM,A)
%     imageschs(DEM,A,pn,pv,...)
%     RGB = imageschs(...)
%
% Description
%
%     Hillshading is a very powerful tool for relief depiction. imageschs
%     calculates a hillshade of a digital elevation model and colors it
%     according to a second grid or matrix. DEM must be an instance of  
%     GRIDobj and A must be a GRIDobj or matrix. 
%
%     imageschs allows to set a number of parameter name/value pairs that
%     control lighting direction and representation of missing values
%     (nan).
%
%     If called with an output variable, imageschs returns an RGB image.
%     The hillshading algorithm follows the logarithmic approach to shaded 
%     relief representation of Katzil and Doytsher (2003).
%
% Input
%
%     DEM         digital elevation model (GRIDobj)
%     A           coloring matrix or GRIDobj
%
% Parameter name/value pairs
%
%     colormap    string for colormap name or [ncol x 3] matrix. Note that 
%                 if NaNs or Infs are found in A, the colormap must not 
%                 have more than 255 colors. Default: 'jet'
%     colorbar    false or true (default)
%     caxis       two element vector defining the value range. Default is
%                 [min(A) max(A)].  
%     truecolor   three element vector (rgb) with values between 0 and 1  
%                 that indicates how true values are plotted if A is 
%                 logical.
%                 Default is [0 1 0].
%     falsecolor  three element vector (rgb) with values between 0 and 1  
%                 that indicates how false values are plotted if A is 
%                 logical.
%                 Default is [1 1 1].
%     nancolor    three element vector (rgb) with values between 0 and 1  
%                 that indicates how NaNs and Infs are plotted 
%                 Default is [1 1 1].
%     azimuth     azimuth angle of illumination, (default=315)
%     altitude    altitude angle of illumination, (default=60)
%     exaggerate  elevation exaggeration (default=2). Increase to
%                 pronounce elevation differences in flat terrain
%     ticklabels  'default', 'nice' or 'none'
%     gridmarkers two element vector with [dx dy] spacing of + markers
%     gridmarkercolor   three element vector (rgb) or color abbreviations
%                 as given in LineSpec (default = 'k')
%                 
%
% Output
%
%     RGB         [DEM.size 3] image (UINT8) with values between 0 and 255
%
% Example
%
%     load exampleDEM
%     H = hillshade(X,Y,dem,315,40);
%     imagesc(X(1,:),Y(:,2),H); axis image; axis xy
%     colormap('gray')
% 
% References
%
%     Katzil, Y., Doytsher, Y. (2003): A logarithmic and sub-pixel approach
%     to shaded relief representation. Computers & Geosciences, 29,
%     1137-1142.
%
% See also: HILLSHADE
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. January, 2013

narginchk(1,inf);
nargoutchk(0,1);

% if A is not supplied to the function, coloring will be according to
% values in DEM
if nargin == 1 || (nargin==2 && isempty(A));
    A = DEM;
end

% Parse inputs
p = inputParser;
p.FunctionName = 'GRIDobj/imageschs';
% required
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'A',@(x) isa(x,'GRIDobj') || ismatrix(A));
% optional
addParamValue(p,'colormap','jet',@(x)(ischar(x) || size(x,2)==3));
addParamValue(p,'caxis',[],@(x) numel(x) == 2);
addParamValue(p,'truecolor',[0 1 0],@(x) isequal(size(x),[1 3]) && (max(x)<=1) && (min(x) >= 0));
addParamValue(p,'falsecolor',[1 1 1],@(x) isequal(size(x),[1 3]) && (max(x(:))<=1) && (min(x(:)) >= 0));
addParamValue(p,'nancolor',[1 1 1],@(x) isequal(size(x),[1 3]) && (max(x(:))<=1) && (min(x(:)) >= 0));
addParamValue(p,'exaggerate',1,@(x) isscalar(x) && x>0);
addParamValue(p,'azimuth',315,@(x) isscalar(x) && x>0);
addParamValue(p,'altitude',60,@(x) isscalar(x) && x>0);
addParamValue(p,'colorbar',true,@(x) isscalar(x));
addParamValue(p,'ticklabels','default',@(x) ischar(x));
addParamValue(p,'gridmarkers',[],@(x) numel(x) == 1 || numel(x) == 2);
addParamValue(p,'gridmarkercolor','k');
parse(p,DEM,A,varargin{:});

% required
DEM        = p.Results.DEM;
A          = p.Results.A;
colmapfun  = p.Results.colormap;
cbar       = p.Results.colorbar;
truecol    = p.Results.truecolor;
falsecol   = p.Results.falsecolor;
exag       = p.Results.exaggerate;
azi        = p.Results.azimuth;
alti       = p.Results.altitude;
nancolor   = p.Results.nancolor;
ticklabels = validatestring(p.Results.ticklabels,{'default','none','nice'});
gridmarkers= p.Results.gridmarkers;
gridmarkercolor = p.Results.gridmarkercolor;

% check if input matrices align
validatealignment(DEM,A)

if isa(A,'GRIDobj')
    A = A.Z;
end

% constrain color range to values given in caxis
if ~isempty(p.Results.caxis)
    A(A<p.Results.caxis(1)) = p.Results.caxis(1);
    A(A>p.Results.caxis(2)) = p.Results.caxis(2);
end

% coordinate matrices
[x,y] = refmat2XY(DEM.refmat,DEM.size);

% nr of colors
nhs = 256;

% calculate hillshading
H = hillshade(DEM,'exaggerate',exag,'azimuth',azi,'altitude',alti);
H = H.Z;
Inan = isnan(H);
if any(Inan(:))
    H(Inan) = 1;
    clear Inan
else
    clear Inan
end
H = gray2ind(H,nhs);

% derive coloring
if ~isa(A,'logical');
%     A(Inan) = nan;
    Inan = isnan(A(:)) | isinf(A(:));
    if any(Inan(:));
        nans = true;
        A(Inan) = nan;
    else
        nans = false;
        clear Inan
    end
    
    if isa(colmapfun,'char');
        ncolors = 256-nans;  
        colmapfun = str2func(lower(colmapfun));
        cmap = colmapfun(ncolors);
    else
        ncolors = size(colmapfun,1);
        if nans && ncolors >= 256;
            error('TopoToolbox:GRIDobj',['NaNs found in the second input argument matrix. \n'...
                  'Please provide colormap with less than 256 colors']);
        else
            cmap = colmapfun;
        end        
    end
    
    if cbar && isempty(p.Results.caxis);
        alims = [min(A(:)) max(A(:))];
    elseif cbar && ~isempty(p.Results.caxis);
        alims = sort(p.Results.caxis,'ascend');
        
    end
    A = gray2ind(mat2gray(A),ncolors);
    
else
    ncolors = 2;
    cmap = [falsecol; truecol];
	nans = false;
    alims = [0 1];
end
    
% create colormap for indexing
cmap = cmap(:);
cmap = bsxfun(@times,cmap,linspace(0,1,nhs));
cmap = reshape(cmap,[ncolors 3 nhs]);
cmap = permute(cmap,[3 1 2]);
cmap = reshape(cmap,[ncolors*nhs 3]);

% create image that indexes into the new colormap
IND  = uint16(H+1) + nhs*uint16(A) + 1;

% handle NaNs
if nans;
    cmapnan   = bsxfun(@times,nancolor,linspace(0,1,nhs)');
    IND(Inan) = uint16(H(Inan)) + nhs*(ncolors) +1;% unclear if this is ok...
    cmap      = [cmap;cmapnan];
end

% same as ind2rgb but returns a mxnx3 matrix with uint8 data
cmapUINT8 = uint8(round(cmap*256));
% see Rob Campbell's mat2im
% http://www.mathworks.de/matlabcentral/fileexchange/26322-mat2im
RGB=reshape(cmapUINT8(IND(:),:),[size(IND),3]);

% plot
if nargout == 0;
    imagesc(x,y,RGB);
    axis xy
    axis image
    
    % add colorbar if needed
    if cbar
        if alims(1) ~= alims(2) ;
            caxis(alims);
        end
        colormap(cmap(nhs:nhs:nhs*ncolors,:));
        colorbar;%('location','south');
    end
    
    % plot nice ticklabels if wanted
    switch ticklabels
        case 'none';
            set(gca,'XTickLabel',{},'YTickLabel',{});
        case 'nice';
            xticklocs = get(gca,'XTick');
            yticklocs = get(gca,'YTick');
            
            set(gca,'XTick',xticklocs([1 end]))
            set(gca,'YTick',yticklocs([1 end]))
            
            set(gca,'XTickLabel',num2str(xticklocs([1 end])','%d'));
            set(gca,'YTickLabel',num2str(yticklocs([1 end])','%d'));
                       
    end
    
    % plot grid
    if ~isempty(gridmarkers);
        if numel(gridmarkers) == 1
            gridmarkers = [gridmarkers gridmarkers];
        end
          
        xgridmarkers = unique(x-rem(x,gridmarkers(1)));
        ygridmarkers = unique(y-rem(y,gridmarkers(2)));
        hold on
        [xx,yy] = meshgrid(xgridmarkers,ygridmarkers);
        plot(xx(:),yy(:),'+','Color',gridmarkercolor);
        hold off
    end
        
        
            
    
elseif nargout == 1;
    rgb = RGB;
end
