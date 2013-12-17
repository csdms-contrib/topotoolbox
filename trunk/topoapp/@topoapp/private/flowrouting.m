function app = flowrouting(hObject,eventdata,app)
% FLOWROUTING creates FLOWobj and STREAMobj from DEM for use in TOPOAPP.
% FLOWROUTING uses by default the 'carving' method to route across sinks.
% For more specific hydrological correction of the digital elevation model 
% (DEM), produce a hydrologically correct DEM first, or modify this code.
% 
% See also: FLOWobj, STREAMobj, imposemin, flowacc, 
% STREAMobj/removeshortstreams
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013

if strcmp(eventdata,'init') % initialize tool
    
    fdicon = imread('fdicon.png','png');
    
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',fdicon,...
        'TooltipString','Set parameters',...
        'Separator','on',...
        'ClickedCallback',{@flowrouting,app});
    
    % add class
    [app] = initclass(app,'STREAMobj','r');
    
else
    
    prompt = {'Minimum upstream area of drainage network [m^2]',...
        'Maximum length of first order streams to remove [m]'};
    answer = inputdlg(prompt,'Set parameters',1,...
        {num2str(app.paras.flowset.minarea),num2str(app.paras.flowset.shorties)});
    
    if ~isempty(answer)
        app.paras.flowset.minarea = str2double(answer{1});
        app.paras.flowset.shorties = str2double(answer{2});
        
        try
            delete(app.gui.Sobj); 
        end
        
        tic
        fprintf(1,'Creating flow direction object...');
        app.FD = FLOWobj(app.DEM,'preprocess','carve');
        app.DEMc = imposemin(app.FD,app.DEM);
        app.FA = flowacc(app.FD);
        fprintf(1,'done. '); toc
        
        fprintf(1,'Creating stream object...');
        app.S = STREAMobj(app.FD,app.FA>(app.paras.flowset.minarea/app.FD.cellsize/app.FD.cellsize));
        if app.paras.flowset.shorties>0
            app.S = removeshortstreams(app.S,app.paras.flowset.shorties);
        end
        fprintf(1,'done. '); toc
        
        hold(app.gui.hax,'on');
        h = findobj(app.gui.hax,'DisplayName','Base STREAMobj');
        if ~isempty(h); delete(h); end

        h = plotobject(app,app.S,'k'); drawnow
        set(h,'DisplayName','Base STREAMobj','visible','off');
        
    end
end
end %