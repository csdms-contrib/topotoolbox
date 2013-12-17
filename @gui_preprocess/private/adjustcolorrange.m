function adjustcolorrange(~,~,app)
% adjust the color range to the data range in the axis
% get axis limits
xlim = get(app.ax,'Xlim');
ylim = get(app.ax,'Ylim');

DEMc = crop(app.DEM,xlim,ylim);
minz = min(DEMc);
maxz = max(DEMc);

if ~isempty(DEMc.Z);
    if minz == maxz;
        set(app.ax,'Clim',[minz-1 minz+1]);
    else
        set(app.ax,'Clim',[minz maxz]);
    end
end
end