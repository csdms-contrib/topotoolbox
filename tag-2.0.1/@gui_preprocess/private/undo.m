function undo(~,~,app)
app.DEM = app.DEMundo;
app.DEMundo = [];
updateimage(app);
set(app.menu.menuundo,'Enable','off');
end   