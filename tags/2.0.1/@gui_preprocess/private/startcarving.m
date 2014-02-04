function startcarving(varargin)
% Sink carving callback
% Callback for the pushbutton.
app = varargin{3};  % Get the structure.


switch get(app.bcarve,'Value')
    case 1
        
        app.DEMundo = app.DEM;
        set(app.menu.menuundo,'Enable','on')
        enabletools(app,'off');
        sethelptext(app,'Define a flowpath by setting two points in downstream direction');
                
        % change pointer
        set(app.fh,'Pointer','crosshair')
        set(app.bcarve,'String','Stop');
        % call carving function
        carvemanually(app);
        
    case 0
        % if finished, recalculate sinks
        app.DEMf = fillsinks(app.DEM);
        
        enabletools(app,'on');
        set(app.bcarve,'String','Start');
        % set pointer to arrow again
        set(app.fh,'Pointer','arrow')
        deletelines(app);
        
end
end