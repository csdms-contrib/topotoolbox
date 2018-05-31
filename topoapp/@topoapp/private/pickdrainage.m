function app = pickdrainage(hObject,eventdata,app)
% PICKDRAINAGE allows interactive drainage extraction within topoapp
%
% See also: topoapp/initclass, topoapp/addobject, STREAMobj/snap2stream, 
% WATERSHEDobj, STREAMobj, removeshortstreams
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013

if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    drainageicon = imread('drainageicon.png','png');
    
    % Set-up toolbar button
    app.gui.TB(end+1) = uitoggletool('Parent',app.gui.hTB,...
        'Cdata',drainageicon,'TooltipString','Pick drainage network',...
        'ClickedCallback',{@pickdrainage,app});
    
    % add classes to topoapp
    [app] = initclass(app,'STREAMobj','r');
    [app] = initclass(app,'WATERSHEDobj','w');
    
else % execute tool
    
    if isempty(app.S)
        warning('No STREAMobj found. Use FLOW routing button first')
        set(hObject,'State','off')
    else
        % Toggles on/off tool to obtain upstream drainages with cursor
        switch get(hObject,'State')
            case 'on'
                set(app.gui.TB,'Enable','off');
                htb = findobj('TooltipString','Pick drainage network');
                set(htb,'Enable','on')
                % change pointer appearance
                enterFcn = @(figHandle, currentPoint) set(app.gui.hfig, 'Pointer','crosshair');
                iptSetPointerBehavior(app.gui.himg,enterFcn);
                iptPointerManager(app.gui.hfig);
                set(app.gui.hfig,'WindowButtonDownFcn',{@PressLeftButton,app});

            case 'off'
                % change pointer appearance
                enterFcn = [];
                iptSetPointerBehavior(app.gui.himg,enterFcn);
                iptPointerManager(app.gui.hfig);
                set(app.gui.hfig,'WindowButtonDownFcn',[]);
                set(app.gui.TB,'Enable','on');
        end
    end
    
end

    function PressLeftButton(~,~,app)
        
        busypointer(app,1)
        % get axis limits
        xlim = get(app.gui.hax,'Xlim');
        ylim = get(app.gui.hax,'Ylim');
        % get current point
        p = get(app.gui.hax,'CurrentPoint');
        pixelx = round(p(1,1));
        pixely = round(p(1,2));
        
        % check if current point is located inside axis
        if pixelx < xlim(1) || pixelx > xlim(2) ...
                || pixely < ylim(1) || pixely > ylim(2)
            % do nothing
        else
            
            [x,y] = sub2coord(app.DEM,pixely,pixelx);
            [~,~,outlet,~,~] = snap2stream(app.S,x,y);
            [DB] = WATERSHEDobj(app.FD,outlet);
            [DN] = STREAMobj(app.FD,'outlets',outlet,'minarea',app.paras.flowset.minarea,'unit','mapunits');
            if app.paras.flowset.shorties>0
                DN = removeshortstreams(DN,app.paras.flowset.shorties);
            end
            
            axes(app.gui.hax), hold on
            h = plotobject(app,DB,app.objects.(class(DB)).color);
            set(h,'Tag',class(DB))
            [app] = addobject(app,DB,'handle',h,'visible',true);
            set(h,'DisplayName',app.objects.(class(DB)).names{end})
            
            if isempty(DN.ix)
                warning('Empty STREAMobj. Watershed may be too small.')
            else
                h = plotobject(app,DN,app.objects.(class(DN)).color);
                set(h,'Tag',class(DN))
                [app] = addobject(app,DN,'handle',h,'visible',true);
                set(h,'DisplayName',app.objects.(class(DN)).names{end})
            end
            
        end
        busypointer(app,0)
    end

end %