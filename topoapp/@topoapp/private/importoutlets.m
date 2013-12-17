function app = importoutlets(hObject,eventdata,app)
% IMPORTOUTLETS creates WATERSHEDobj and STREAMobj from outlet
% coordinate pairs (lat,lon) imported to topoapp
%
% See also: snap2stream, WATERSHEDobj, STREAMobj, addobject, initclass
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    outletsicon = imread('outletsicon.png','png');
    
    % Set up toolbar button
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',outletsicon,'TooltipString','Import outlets from ascii file (X,Y)',...
        'ClickedCallback',{@importoutlets,app});
    
    % add classes to topoapp
    [app] = initclass(app,'STREAMobj','r');
    [app] = initclass(app,'WATERSHEDobj','w');
    
else % execute tool
    
    if isempty(app.S)
        error('No STREAMobj found. Use FLOW routing button first')
    end
    
    % Import data
    S = uiimport;
    fieldname = fieldnames(S);
    outlets = S.(fieldname{1});
    
    if max(abs(outlets(:,1)))<=90 && max(abs(outlets(:,2)))<=180 % probably lat,lon
        fprintf(1,'\n   Lat\t\tLon\n');
        disp(outlets)
        fprintf(1,'Outlets appear to be in geographic coordinates.\n');
        fprintf(1,'\nTrying to convert to UTM coordinates.\n');
        try
            [outlets(:,2),outlets(:,1)] = projfwd(app.DEM.georef,outlets(:,1),outlets(:,2));
            fprintf(1,'\n   X\t\tY\n');
            disp(outlets)
        catch
            error('Problems converting lat, lon to projected coordinates. Stopping import.')
        end
    end
    
    % Snap outlets to streams
    fprintf(1,'\nSnapping outlets to STREAMobj...');
    dlg_title = 'Input for snap2stream';
    prompt = {'Maximum distance (map units)',...
        'Stream order (e.g., >3, ==5, etc.)',...
        'Snap to feature (all | confluence | outlet | channelhead)'};
    answer = inputdlg(prompt,dlg_title,1,{'inf','','all'});
    [~,~,outlets_ix,~,~] = snap2stream(app.S,outlets(:,2),outlets(:,1),...
        'snapto',answer{3},'maxdist',str2double(answer{1}),'streamorder',answer{2});
    fprintf(1,'done.\n');
    
    % Watershed identification
    wschoice = questdlg('Find watersheds:','Watershed identification', ...
        'Automatically','Interactively','Cancel','Automatically');
    
    if strcmp(wschoice,'Cancel'); return; end
    
    for i = 1 : length(outlets_ix)
        
        fprintf(1,'\nSearching watershed for X = %1.4f, Y = %1.4f  ',outlets(i,2),outlets(i,1));
        if isnan(outlets_ix(i))
            warning('Outlet appears to be located outside of current map extent.')
        else
            [DB] = WATERSHEDobj(app.FD,outlets_ix(i));
            fprintf(1,'Area = %1.0f m^2\n',numel(DB.ix)*app.DEM.cellsize*app.DEM.cellsize);
            axes(app.gui.hax); hold on
            h = plotobject(app,DB,app.objects.(class(DB)).color);
            
            if strcmp(wschoice,'Interactively')
                
                set(h,'Selected','on')
                reply = input('Modify watershed? Y/N [N]: ', 's');
                if isempty(reply); reply = 'N'; end
                reply = upper(reply);
                set(h,'Selected','off')
                
                while strcmp(reply,'Y')
                    % get axis limits
                    xlim = get(app.gui.hax,'Xlim');
                    ylim = get(app.gui.hax,'Ylim');
                    % get current point
                    axes(app.gui.hax);
                    fprintf(1,'\nPick outlet from DEMfigure. ');
                    [pixelx,pixely] = ginput(1);
                    % check if current point is located inside axis
                    if pixelx < xlim(1) || pixelx > xlim(2) ...
                            || pixely < ylim(1) || pixely > ylim(2)
                        % do nothing
                    else
                        [x,y] = sub2coord(app.DEM,pixely,pixelx);
                        fprintf(1,'X = %1.4f, Y = %1.4f ',x,y);
                        [~,~,outlets_ix(i),~,~] = snap2stream(app.S,x,y);
                        [DB] = WATERSHEDobj(app.FD,outlets_ix(i));
                        fprintf(1,'Area = %1.0f m^2\n',numel(DB.ix)*app.DEM.cellsize*app.DEM.cellsize);
                        axes(app.gui.hax), hold on
                        delete(h)
                        h = plotobject(app,DB,app.objects.(class(DB)).color);
                        set(h,'Selected','on')
                    end
                    reply = input('Modify watershed? Y/N [N]: ', 's');
                    if isempty(reply); reply = 'N'; end
                    reply = upper(reply);
                    set(h,'Selected','off')
                end
            end
            
            set(h,'Tag',class(DB))
            [app] = addobject(app,DB,'handle',h,'visible',true);
            set(h,'DisplayName',app.objects.(class(DB)).names{end})
            [DN] = STREAMobj(app.FD,'outlets',outlets_ix(i),'minarea',app.paras.flowset.minarea,'unit','mapunits');
            if app.paras.flowset.shorties>0
                DN = removeshortstreams(DN,app.paras.flowset.shorties);
            end
            
            if isempty(DN.ix)
                warning('Empty STREAMobj. Watershed may be too small.')
            else
                h = plotobject(app,DN,app.objects.(class(DN)).color);
                set(h,'Tag',class(DN))
                [app] = addobject(app,DN,'handle',h,'visible',true);
                set(h,'DisplayName',app.objects.(class(DB)).names{end})
            end
        end
    end
    hold off
    fprintf(1,'done.\n');
    
    
end




end %




