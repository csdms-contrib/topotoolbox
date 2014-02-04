function app = streamview(hObject,eventdata,app)
% STREAMVIEW toggles displaying of base STREAMobj in topoapp
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    flowicon = imread('flowicon.png','png');    
        
    % Set up toolbar button
    app.gui.TB(end+1) = uitoggletool('Parent',app.gui.hTB,...
        'Cdata',flowicon,...
        'TooltipString','Show Stream network',...
        'ClickedCallback',{@streamview,app});
    
else
    
    if isempty(app.S)
        warning('No STREAMobj found. Use FLOW routing button first')
    else
        switch get(hObject,'State')
            case 'on'
                h = findobj(app.gui.hax,'Displayname','Base STREAMobj');
                set(h,'visible','on')
            case 'off'
                h = findobj(app.gui.hax,'Displayname','Base STREAMobj');
                set(h,'visible','off')
        end
    end

end
end %