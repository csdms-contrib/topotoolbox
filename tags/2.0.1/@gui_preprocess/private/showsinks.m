function showsinks(~,~,app)
% switch between sink display

state = get(app.menu.menuss,'Checked');

switch state
    case 'off'
        set(app.menu.menuss,'Checked','on');        
        set(app.im,'AlphaDataMapping','none','AlphaData',1-app.SINKS.Z*.5);    

    case 'on'
        set(app.menu.menuss,'Checked','off');
        set(app.im,'Cdata',app.DEM.Z,'AlphaData',1); 
        
end
end