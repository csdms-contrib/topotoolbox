function varargout = densify(S,spacing)

%DENSIFY Increase number of vertices in stream network using splines
%
% Syntax
%
%     [x,y] = densify(S,spacing)
%     MS    = densify(S,spacing)
%
% Description
%
%     This function increases the number of vertices in a stream network
%     (STREAMobj) using spline interpolation. The function returns either
%     nan-separated lists of x- and y-coordinates or a mapstruct (MS) which
%     can be exported to a shapefile using the function shapewrite
%     (requires the mapping toolbox).
%
% Input arguments
%
%     S        STREAMobj
%     spacing  spacing of vertices in the resampled network (e.g. 5 [m])
%
% Output arguments
%
%     x,y      nan-separated coordinate vectors
%     MS       mapstruct
%
%
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 16. Octobre, 2014



if nargout == 2
    I = streampoi(S,'co','logical');
    [x,y,d,i] = STREAMobj2XY(S,S.distance,I);
    
    ixe = find(isnan(x)) - 1;
    ixe = ixe(:);
    ixs = [1; ixe(1:end-1)+2];
    
    xn = [];
    yn = [];
    for r = 1:numel(ixs);
        % vector with interpolation locations
        dmax = d(ixs(r));
        dmin = d(ixe(r));
                
        dd   = linspace(dmin,dmax,max(ceil((dmax-dmin)/spacing),2));
        
        % vector with confluence locations
        dc   = d(ixs(r):ixe(r));
        ii   = i(ixs(r):ixe(r));
        dc   = dc(logical(ii));
        
        dd   = unique([dd(:);dc(:)],'sorted');
        
        xxyy = splinedensify(d(ixs(r):ixe(r)),x(ixs(r):ixe(r)),y(ixs(r):ixe(r)),dd);
        
        xn   = [xn; xxyy(:,1); nan];
        yn   = [yn; xxyy(:,2); nan];
        
    end
    varargout{1} = xn;
    varargout{2} = yn;
    
    
else
    
    MS = STREAMobj2mapstruct(S);
    
    for r = 1:numel(MS);
        x = MS(r).X;
        y = MS(r).Y;
        
        x = x(1:end-1);
        y = y(1:end-1);
        
        x = x(:);
        y = y(:);
        d = [0; cumsum(sqrt(diff(x).^2 + diff(y).^2))]';
        
        dd = 0:spacing:d(end);
        
        if dd(end) ~= d(end);
            dd(end+1) = d(end);
        end
        
        xxyy = splinedensify(d(:),x(:),y(:),dd(:));
        MS(r).X = [xxyy(:,1);nan];
        MS(r).Y = [xxyy(:,2);nan];
        
    end
    varargout{1} = MS;
end
        
    
        

end





function xxyy = splinedensify(d,x,y,dd)

xxyy = spline(d,[x y]',dd);
xxyy = xxyy';
end
