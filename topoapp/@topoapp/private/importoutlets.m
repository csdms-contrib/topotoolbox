function app = importoutlets(hObject,eventdata,app)
% IMPORTOUTLETS creates WATERSHEDobj and STREAMobj from outlet
% coordinates imported to topoapp
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
        warning('No STREAMobj found. Use FLOW routing button first')
    else
        [filename,pathname] = uigetfile('*.txt','Import outlet data');

        if filename
            outlets = importdata([pathname,filename]);
            if max(outlets(:,1))<=180 && max(outlets(:,2))<=90 % probably lat,lon
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
            axes(app.gui.hax); hold on
            fprintf(1,'\nSearching for watersheds...');
            [~,~,outlets,~,~] = snap2stream(app.S,outlets(:,2),outlets(:,1));
            for i = 1 : length(outlets)

                if isnan(outlets(i))
                    warning('Outlet appears to be located outside of current map extent.')
                else
                    [DB] = WATERSHEDobj(app.FD,outlets(i));
                    [DN] = STREAMobj(app.FD,'outlets',outlets(i),'minarea',app.paras.flowset.minarea,'unit','mapunits');
                    if app.paras.flowset.shorties>0
                        DN = removeshortstreams(DN,app.paras.flowset.shorties);
                    end

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
                        set(h,'DisplayName',app.objects.(class(DB)).names{end})
                    end
                end
            end
            hold off
            fprintf(1,'done.\n');
        end
    end
end

end %




