function startfilling(varargin)

% Sink filling callback
% Callback for the pushbutton.
app = varargin{3};  % Get the structure.


switch get(app.bfill,'Value')
    
    % button is on
    case 1
        set(app.bfill,'String','Stop');
        app.DEMundo = app.DEM;
        set(app.menu.menuundo,'Enable','on')
        
        % which option is chosen
        meth = find(cell2mat(get(app.fillmeth,'Value')));
        % if method is 1 or 2;
        if ismember(meth,[1 2]);
            % set helptext
            sethelptext(app,'Fill all sinks with predefined conditions.')
            
            % get value from edit text box
            val = str2double(get(app.fillvaledit,'string'));
            if meth == 1;
                app.DEMf = fillsinks(app.DEM);
                I    = (app.DEMf.Z-app.DEM.Z) > 0;
                I    = I & ~bwareaopen(I,val+1);
                app.DEM.Z(I) = app.DEMf.Z(I);
            elseif meth == 2;
                DEMf = fillsinks(app.DEM,val);
                app.DEM = DEMf;                
                
            end
                      
            updateimage(app) 
            sethelptext(app,'Sinks are filled.')           
            set(app.bfill,'Value',0);
            set(app.bfill,'String','Start');
           
        elseif meth == 3;
            % make tools unabled
            app.DEMf = fillsinks(app.DEM);
            enabletools(app,'off')
            sethelptext(app,'Click on sinks that should be filled.')
            set(app.fh,'Pointer','crosshair')            
            fillmanually(app);
            
        end
        
        
    case 0
        % turn off filling procedure
        set(app.bfill,'String','Start');
        enabletools(app,'on')
        set(app.fh,'Pointer','arrow')
end
end