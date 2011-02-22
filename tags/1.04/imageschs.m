function varargout = imageschs(X,Y,dem,A,exagg,B)

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

h = gca; %('DataAspectRatio',[1 1 1]);

g = image('Parent',h,...
          'XData',X(1,:),'YData',Y(:,1),'Cdata',+A,...
          'CdataMapping','scaled',...
          'AlphaData',H);
    
axis image;
axis xy;

if ~ishold(h)
    hold off
end

% g = image(X(1,:),Y(:,1),+A);
% axis image;
% axis xy;

% set(g,'AlphaData',H);



if nargin == 6;
    x = X(1,:);
    y = Y(:,1);
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(x(boundary(:,2)), y(boundary(:,1)), 'w', 'LineWidth', 2)
    end
    hold off
end

if nargout == 1;
    varargout{1} = h;
end


