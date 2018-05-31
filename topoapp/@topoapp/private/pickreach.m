function app = pickreach(hObject,eventdata,app)
% PICKREACH allows interactive extraction of a reach within topoapp
%
% See also: topoapp/initclass, topoapp/addobject, 
% STREAMobj/STREAMobj2GRIDobj, REACHobj
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013

if strcmp(eventdata,'init') % initialize tool
    
    % Load button icons
    reachicon = imread('reachicon.png','png');
    
    % Set up toolbar button
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',reachicon,'TooltipString','Pick channel reach',...
        'ClickedCallback',{@pickreach,app});
    
    % add class to topoapp
    app = initclass(app,'REACHobj','g');
    
else % execute tool
    
    if isempty(app.S)
        warning('No STREAMobj found. Use FLOW routing button first')
    else
        % Define reach along channels with two points
        set(app.gui.TB,'Enable','off');
        W = STREAMobj2GRIDobj(app.S);
        [x,y] = ginput(1); % Let user pick first point
        ix1 = sub2ind(app.FD.size,round(y(1)),round(x(1)));
        hold on; ht(1) = plot(x,y,'ro'); hold off
        [x,y] = ginput(1); % Let user pick second point
        ix2 = sub2ind(app.FD.size,round(y(1)),round(x(1)));
        hold on; ht(2) = plot(x,y,'ro'); hold off

        [ix1,~] = snap2stream(W,ix1);
        [ix2,~] = snap2stream(W,ix2);

        R = REACHobj(app.FD,[ix1,ix2]);

        axes(app.gui.hax), hold on
        h = plotobject(app,R,app.objects.(class(R)).color); hold off
        set(h,'Tag',class(R))
        [app] = addobject(app,R,'handle',h,'visible',true);
        set(h,'DisplayName',app.objects.(class(R)).names{end})
        set(app.gui.TB,'Enable','on');
        delete(ht)
        figure, plotzobj(app,R,false,0);
        title(app.objects.(class(R)).names{end})
        set(app.gui.TB,'Enable','on');
    end
    
end

end %