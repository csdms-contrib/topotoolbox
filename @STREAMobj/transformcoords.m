function varargout = transformcoords(S,varargin)

% transform coordinates of stream network
%
% Syntax
%
%     St = transformcoords(S)
%     St = transformcoords(S,'pn',pv,...)
%     [xt,yt] = ...
%
% Description
%
%     transformcoords lets you transform the spatial coordinates. By default
%     the transformation will be done interactively by drawing a new x-axis 
%     into a planform layout of the river network. Alternatively, the
%     coordinate transformation will be performed automatically using a
%     principal component analysis.
%
% Input arguments
%     
%     S      STREAMobj
%     
%     Parameter name/value pairs {default}
%
%     'method'         {'interactive'} or 'pca'
%     'nonnegcoords'   true or {false}. Setting to true will set minimum x
%                      and y values to zero.
%
% Output arguments
%
%     St     transformed STREAMobj
%     xt,yt  node-attribute lists of transformed coordinates
%
%
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 23. May, 2017



p = inputParser;
p.FunctionName = 'STREAMobj/transformcoords';
addParameter(p,'method','interactive');
addParameter(p,'nonnegcoords',false);
parse(p,varargin{:});

method = validatestring(p.Results.method,{'pca','interactive','inv'});

switch method
    case 'pca'
        [~,s] = pca([S.x S.y]);
        xx = s(:,1);
        yy = s(:,2);
    case 'interactive'
        plot(S);
        hold on
        h = imline;
        title('Double-click line to finish')
        position = wait(h);
        
        sp = position(1,:);
        ep = position(2,:);
        
        x  = S.x-sp(1);
        y  = S.y-sp(2);
        
        ep = ep-sp;
        ep = ep(:);
        
        % projection of [x; y] on vector a
        beta = (ep' * [x y]')./(ep'*ep);
        xx = beta .* sqrt(ep'*ep);
        yy = beta .* ep;
        yy = [x y]-yy';
        
        s  = sign(yy(:,2));
        yy = hypot(yy(:,1),yy(:,2));
        yy = yy.*s;
    case 'inv'
        nrrows = S.size(1);
        nrcols = S.size(2);

        x = [ones(nrcols,1) (1:nrcols)' ones(nrcols,1)]*S.refmat;
        xx = x(:,1)';

        y = [(1:nrrows)' ones(nrrows,2)]*S.refmat;
        yy = y(:,2)';
end

if p.Results.nonnegcoords
    xx = xx - min(xx);
    yy = yy - min(yy);
end

if nargout == 1
    S.x = xx(:);
    S.y = yy(:);
    varargout{1} = S;
    
else
    varargout{1} = xx;
    varargout{2} = yy;
    
end
        


