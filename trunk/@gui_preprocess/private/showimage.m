function showimage(~,~,app)
set(app.menu.menuim,'Checked','on')
set(app.menu.menuhs,'Checked','off')
set(app.im,'CData',app.DEM.Z)
end