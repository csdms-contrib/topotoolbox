function varargout = imageschs(X,Y,dem,A,exagg,boundaries)

% display image with hillshading
%
% Syntax
%
%     imageschs(X,Y,dem,A)
%     imageschs(X,Y,dem,A,exagg)
%     h = ...
%
% Description
%
%     imageschs is a wrapper for the function imagesc. In addition, it
%     applies a hillshading to the display of the matrix/image A based on
%     the digital elevation model dem with the coordinates created by
%     meshgrid X,Y. X, Y, dem and A must have same size. exagg (by default
%     = 1) is the vertical exaggeration of the DEM.
%     h = imageschs(...) returns the handle for an image graphics object.
%
% Example
%
%     load exampleDEM
%     A = ezflowacc(X,Y,dem);
%     imageschs(X,Y,dem,A,3)
%
% See also: HILLSHADE, IMAGESC
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 18. May, 2009


error(nargchk(3, 6, nargin))



if isvector(X) && isvector(Y)
    [X,Y] = meshgrid(X,Y);
end

% check if Y is increasing
if Y(1)>Y(2);
    Y = flipud(Y);
    dem = flipud(dem);
    if nargin > 3;
    A   = flipud(A);
    end
end

if nargin == 3;
    hillshade(X,Y,dem); 
    if nargout>0;
    varargout{1} = gca;
    end
    return
elseif nargin == 4;
    exagg = 1;
end
    
if ~isequal(size(X),size(Y),size(dem),size(A))
    error('all matrices must have same size')
end



% calculate hillshading
H = hillshade(X,Y,dem,[],[],exagg);

% check axis

clf('reset')
h = gca; %('DataAspectRatio',[1 1 1]);
cla(h,'reset')


I       = isnan(A) | isinf(A);
naninf  = any(I(:));
if naninf
    indexed = true;
else
    indexed = false;
end

if indexed;
    A_min   = min(A(~I(:)));
    A_max   = max(A(~I(:)));
    A_range = A_max-A_min;
    if A_range == 0;
        A_range = 1;
    end
    
    nrA     = 256 - double(naninf);
    B       = (A-A_min)/A_range;
    B       = B*nrA + naninf;
    
    if naninf
        B(I) = 0;
        B       = uint8(B);
        colormap([0 0 0; jet(nrA)]);
    else
        B       = uint8(B);
        clear I
        colormap(jet(nrA))
    end
        
    
    g = image('Parent',h,...
          'XData',X(1,:),'YData',Y(:,1),'Cdata',B,...
          'CdataMapping','direct',...
          'AlphaData',H);
    cbar_axes = colorbar;
    axis(cbar_axes,[0 1 1 nrA]);
    cbar_labels = get(cbar_axes,'Yticklabel');
    cbar_labels = str2num(cbar_labels); %#ok<ST2NM>
    cbar_labels = (cbar_labels/(nrA - double(naninf)) * A_range)+A_min;
    set(cbar_axes,'Yticklabel',cbar_labels);
    
    box on;
    
    
else
    g = image('Parent',h,...
          'XData',X(1,:),'YData',Y(:,1),'Cdata',+A,...
          'CdataMapping','scaled',...
          'AlphaData',H,...
          'Clipping','on');
    box on
    
    if islogical(A);
        colormap([0 0 0; 1 0 0]);
    end
end
    
axis image;
axis xy;

if ~ishold(h)
    hold off
end

if nargin == 6;
    x = X(1,:);
    y = flipud(Y(:,1));
    hold on
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        plot(x(boundary(:,2)), y(boundary(:,1)), 'k', 'LineWidth', 1)
    end
    hold off
end

if nargout == 1;
    varargout{1} = h;
end


