function lcpcarving(app,IX)
% Least cost carving
% IX contains are two linear indices into the DEM
% IX(1) = start index
% IX(2) = end index

val = str2double(get(app.carvevaledit,'string'));
% create a mask. Within the mask the least cost path will be found
[xy(:,1) xy(:,2)]  = ind2coord(app.DEM,IX);
xyse = xy;

% find radius center
c = sum(xy,1)/2;
% radius
dxy = xy(1,:)-xy(2,:);
r = hypot(dxy(1),dxy(2))/2;

% bounding box
xy(1,:) = c + val*[r r];
xy(2,:) = c - val*[r r];

% take care that bounding box doesn't extend beyond dem limits
xy(:,1) = max(min(xy(:,1),max(app.X)),min(app.X));
xy(:,2) = max(min(xy(:,2),max(app.Y)),min(app.Y));

% coordinate to index
IX  = coord2ind(app.DEM,...
               [min(xy(:,1)) max(xy(:,1))],[min(xy(:,2)) max(xy(:,2))]);
           
hold(app.ax,'on')
g = plot(xy([1 2 2 1 1],1),xy([1 1 2 2 1],2),'-k');
pause(.1)
hold(app.ax,'off')

Z = crop(app.DEM,xy(:,1),xy(:,2));
IXSE = coord2ind(Z,xyse(:,1),xyse(:,2));

DD  = graydist((Z.Z-min(Z.Z(:))+1).^2,IXSE(2),'quasi-euclidean');
D   = Z;
D.Z = DD;

FD     = FLOWobj(D);
IXpath = flowpathextract(FD,IXSE(1));

for r = 2:numel(IXpath)
    Z.Z(IXpath(r)) = min(Z.Z(IXpath(r-1)),Z.Z(IXpath(r)));
end
        
pause(.4);
delete(g);

writeinmatrix(app,Z,IX);

[xx,yy] = ind2coord(Z,IXpath);
hold(app.ax,'on')
app.lh(end+1) = plot(xx,yy,'k-');
pause(.1)
hold(app.ax,'off')

end