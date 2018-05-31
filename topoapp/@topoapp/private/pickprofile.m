function app = pickprofile(hObject,eventdata,app)
% PICKPROFILE allows interactive extraction of a profile within topoapp. 
% Double-klick ends point selection
%
% See also: topoapp/initclass, topoapp/addobject, getline
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013

if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    profileicon = imread('profileicon.png','png');
    
    % Set up toolbar button
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',profileicon,'TooltipString','Draw profile',...
        'ClickedCallback',{@pickprofile,app});
    
    % add class to topoapp
    [app] = initclass(app,'PROFILEobj','k');
    
else % execute tool
    
    axes(app.gui.hax)
    set(app.gui.TB,'Enable','off');
    [x,y] = getline(app.gui.hax);
    
    P = PROFILEobj(app.DEM,x,y);
    hold on
    h = plotobject(app,P,app.objects.(class(P)).color); hold off
    set(h,'Tag',class(P))
    [app] = addobject(app,P,'handle',h,'visible',true);
    set(h,'DisplayName',app.objects.(class(P)).names{end})
    set(app.gui.TB,'Enable','on');
    
    figure, plotzobj(app,P,true,0);
    
end

end %





