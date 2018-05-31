function h = plotobject(app,objdata,color)
% PLOTOBJECT creates simple map-view plots of topoapp objects
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013
if isa(objdata,'STREAMobj')
    nrc = numel(objdata.x);
    M = sparse(double(objdata.ix),double(objdata.ixc),true,nrc,nrc);
    [R,C] = ind2sub(app.DEM.size,double(objdata.IXgrid));
    [x,y] = gplot(M,[C R]);
    h = plot(app.gui.hax,x,y,'color',color);
    
elseif isa(objdata,'SWATHobj')
    SW = objdata;
    for i = 1 : length(SW.X)
        SW.X{i} = interp1(app.X,(1:app.DEM.size(2)),SW.X{i});
        SW.Y{i} = interp1(app.Y,(1:app.DEM.size(1)),SW.Y{i});
    end
    x = cell2mat(SW.X);
    y = cell2mat(SW.Y);
    x = [x; nan(1,size(x,2))];
    y = [y; nan(1,size(y,2))];
    x = reshape(x,numel(x),1);
    y = reshape(y,numel(y),1);
    h = plot(app.gui.hax,x,y,[color,'o'],'Markersize',2);
    
elseif not(isempty(objdata))
    [y,x] = coord2sub(app.DEM,objdata.x,objdata.y);
    h = plot(app.gui.hax,x,y,'color',color);
    
else
    % do nothing
    h = [];
end
end %