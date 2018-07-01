function app = demview(hObject,eventdata,app)
% DEMVIEW toggles between different DEM displays in TOPOAPP
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013

if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    demicon = imread('demicon.png','png');   
        
    % Set up toolbar button
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',demicon,...
        'TooltipString','Toggle DEM display',...
        'Userdata',2,...
        'ClickedCallback',{@demview,app});
    
else
    
    switch get(hObject,'Userdata')
        
        case 1
            % display grey-scale hillshade
            set(app.gui.himg,'Cdata',app.HS,'Cdatamapping','scaled');
            colormap(gray);
            set(hObject,'Userdata',2)
        case 2
            % display elevation-colored hillshade
            set(app.gui.himg,'Cdata',app.RGB)
            colormap(jet)
            set(hObject,'Userdata',3)
            
        case 3
            % display slope-colored hillshade
            set(app.gui.himg,'Cdata',app.HSG)
            colormap(jet)
            set(hObject,'Userdata',4)
            
        case 4
            % display sink depths
            set(app.gui.himg,'Cdata',app.SINKS)
            colormap(jet)
            set(hObject,'Userdata',1)
    end

end
end %