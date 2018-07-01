function app = swathmapping(hObject,eventdata,app)
% SWATHMAPPING allows creating SWATHobj from topoapp objects in topoapp
%
% See also: topoapp/listfigure, topoapp/initclass, topoapp/addobject,
% SWATHobj, SWATHobj/mapswath, SWATHobj/modify, SWATHobj/plot,
% SWATHobj/plotdz, SWATHobj/tidyswath
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    swathicon = imread('swathicon.png','png');
    
    % Set-up toolbar button
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',swathicon,'TooltipString','Swath mapping',...
        'ClickedCallback',{@swathmapping,app});
    
    % add class
    [app] = initclass(app,'SWATHobj','m');
    
else % execute tool
    
    SW = [];
    paras = defaultparas;
    
    % Initialize basic figure
    definedclasses = {'PROFILEobj','STREAMobj','WATERSHEDobj','REACHobj'};
    [hf,posYlast] = listfigure(app,definedclasses);
    set(findobj('Tag','Itemlist'),'Max',1); %  limit to one object at a time
    % Add tool-specific options
    addoptions(hf,posYlast,app);
    
end

    %------------------------------------------------------------
    function addoptions(hfig,posYlast,app)
        
        figure(hfig)
        %--------------------------------------------------------
        % positioning:
        dLs = 0.025; dLl = 0.2833;
        Lft1 = dLs; Lft2 = Lft1+(dLl+dLs+dLs); Lft3 = Lft2+(dLl+dLs+dLs);
        dLftedit1 = 0.12;  dLftedit2 = 0.2133;
        Hpb = 0.06; Hedit = 0.045; Htxt = 0.045;
        Wpb = dLl; Wedit = 0.07; Wtxts = 0.1; Wtxtl = 0.17;
        %--------------------------------------------------------
        
        %% Parameter controls
        % Link to streamobj
        posY = posYlast; ct = 1;
        h(ct) = uicontrol('Style','checkbox',...
            'units','normalized','position',[Lft3,posY,Wpb,Htxt],...
            'String','Link distance to STREAMobj',...
            'TooltipString','Only available for reaches',...
            'Callback',{@link2streamobj,app});
        posY = posY - 0.05; ct = ct+1;
        % Swath width
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','Position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String',num2str(paras.width)); ct = ct+1;
        h(ct) = uicontrol('Style','text',...
            'String','Swath width (m)',...
            'units','normalized','Position',[Lft3,posY,Wtxtl,Htxt],...
            'Userdata',h(ct-1)); ct = ct+1;
        posY = posY - 0.05;
        % Gap width
        h(ct) = uicontrol('Style','edit',...
            'String',num2str(paras.gap),...
            'units','normalized','Position',[Lft3+dLftedit2,posY,Wedit,Hedit]); ct = ct+1;
        h(ct) = uicontrol('Style','text',...
            'String','Central gap (m)',...
            'units','normalized','Position',[Lft3,posY,Wtxtl,Htxt],...
            'Userdata',h(ct-1)); ct = ct+1;
        posY = posY - 0.05;
        % Left/right
        h(ct) = uicontrol('Style','checkbox',...
            'String','Left','Value',1,...
            'units','normalized','Position',[Lft3,posY,2*Wedit,Hedit]); ct = ct+1;
        h(ct) = uicontrol('Style','checkbox',...
            'String','Right','Value',1,...
            'units','normalized','Position',[Lft3+dLl/2,posY,2*Wedit,Hedit]); ct = ct+1;
        posY = posY - 0.05;
        % Resample dx checkbox
        h(ct) = uicontrol('Style','edit',...
            'String',num2str(paras.stepx),...
            'units','normalized','Position',[Lft3+dLftedit2,posY,Wedit,Hedit]); ct = ct+1;
        h(ct) = uicontrol('Style','checkbox','value',1,...
            'Tag','swaths.resamplex',...
            'String','Resample x (m)',...
            'units','normalized','Position',[Lft3,posY,Wtxtl,Hedit],...
            'Userdata',h(ct-1)); ct = ct+1;
        posY = posY - 0.05;
        % Resample dy checkbox
        h(ct) = uicontrol('Style','edit',...
            'String',num2str(paras.stepx),...
            'units','normalized','Position',[Lft3+dLftedit2,posY,Wedit,Hedit]); ct = ct+1;
        h(ct) = uicontrol('Style','checkbox','value',1,...
            'String','Resample y (m)',...
            'units','normalized','Position',[Lft3,posY,Wtxtl,Hedit],...
            'Userdata',h(ct-1)); ct = ct+1;
        posY = posY - 0.05;
        % Smoothing (slider max = max possible for filtfilt)
        h(ct) = uicontrol('Style','slider',...
            'units','normalized','Position',[Lft3,posY-0.05,Wpb,Hedit],...
            'Min',1,'Max',paras.maxsmooth,...
            'Value',paras.smoothvalue); ct = ct+1;
        h(ct) = uicontrol('Style','checkbox',...
            'String','Smoothing (set below)',...
            'units','normalized','Position',[Lft3,posY,Wpb,Hedit],...
            'Userdata',h(ct-1)); ct = ct+1;
        posY = posY - 0.1;
        % Keep nodes
        h(ct) = uicontrol('Style','checkbox',...
            'String','Keep nodes',...
            'units','normalized','Position',[Lft3,posY,Wpb,Hedit],...
            'Value',0); ct = ct+1;
        posY = posY - 0.05;
        % Keep distance
        h(ct) = uicontrol('Style','checkbox',...
            'String','Keep distance',...
            'units','normalized','Position',[Lft3,posY,Wpb,Hedit],...
            'Value',0);
        posY = posY - 0.05; ct = ct+1;
        % Remove duplicates
        h(ct) = uicontrol('Style','checkbox',...
            'String','Remove duplicates',...
            'units','normalized','Position',[Lft3,posY,Wpb,Hedit],...
            'Value',0); ct = ct+1;
        posY = posY - 0.05;
        % Show outline
        h(ct) = uicontrol('Style','checkbox',...
            'String','Show outline',...
            'units','normalized','Position',[Lft3,posY,Wpb,Hedit],...
            'Value',1); ct = ct+1;
        posY = posY - 0.05;
        % Show points
        h(ct) = uicontrol('Style','checkbox',...
            'String','Show points',...
            'units','normalized','Position',[Lft3,posY,Wpb,Hedit],...
            'Value',0);
        paract = ct; ct = ct+1;
        %--------------------------------------------------------
        
        %% Tool options
        % Plot swath elevation statistics
        h(ct) = uicontrol('Style','pushbutton',...
            'String','Map swath',...
            'units','normalized','Position',[Lft2,.275,Wpb,Hpb],...
            'Callback',{@mkswath,app},...
            'Userdata',h(1:paract)); ct = ct+1;
        h(ct) = uicontrol('Style','pushbutton',...
            'String','Plot elevation swath',...
            'units','normalized','Position',[Lft2,.2,Wpb,Hpb],...
            'Callback',{@plotswath,app},...
            'Userdata',h(1:paract)); ct = ct+1;
        % Plot swath slope statistics
        h(ct) = uicontrol('Style','pushbutton',...
            'String','Plot slope swath',...
            'units','normalized','Position',[Lft2,.125,Wpb,Hpb],...
            'Callback',{@plotswath,app},...
            'Userdata',h(1:paract)); ct = ct+1;
        % Add swath to topoapp
        h(ct) = uicontrol('Style','pushbutton',...
            'String','Add to topoapp',...
            'units','normalized','Position',[Lft2,.05,Wpb,Hpb],...
            'Callback',{@addswath,app},...
            'Userdata',h(1:paract));
        set(h,'Interruptible','off'); clear h;
        set(findobj('Style','text'),'Backgroundcolor',[.8 .8 .8])
        set(findobj('Style','checkbox'),'Backgroundcolor',[.8 .8 .8])
        %--------------------------------------------------------
        
    end
    %------------------------------------------------------------
    function paras = defaultparas
        
        % swath-mapping parameters
        paras.width = 2e4;
        paras.gap = 0;
        paras.stepx = 500;
        paras.stepy = 500;
        paras.smooth = 0;
        paras.smoothvalue = 10;
        paras.maxsmooth = 10000;
        paras.showoutline = 1;
        paras.showpoints = 0;
        paras.keepnodes = 0;
        paras.keepdistance = 0;
        %paras.curvthreshold = 0.002;
        
        paras.mincol = 0;
        paras.maxcol = 60;
        paras.valleybot = 3;

    end % 
    %------------------------------------------------------------
    function link2streamobj(~,~,app)
        if get(gcbo,'Value')
            if ~isempty(app.objects.STREAMobj.data)
                % select streamobj
                strobjname = app.objects.STREAMobj.names;
                set(0,'DefaultUIcontrolUnits','pixel');
                [s,ok] = listdlg('PromptString','Select STREAM object to link with:',...
                                'SelectionMode','single',...
                                'ListString',strobjname);
                set(0,'DefaultUIcontrolUnits','normalized');
                if ok
                    set(gcbo,'Userdata',s);
                else
                    set(gcbo,'Value',0)
                end
            else
                fprintf(1,'No STREAMobj to link with.\n');
                set(gcbo,'Value',0)
            end
        else
            set(gcbo,'Userdata',[]);
        end
    end % 
    %------------------------------------------------------------
    function [offset] = adjustdistance(classname,ix,app)
        offset = 0;            
        h = findobj('String','Link distance to STREAMobj','Visible','on');
        if get(h,'Value') & ( strcmp(classname,'STREAMobj') || ...
            strcmp(classname,'REACHobj') )
        
            iy = get(h,'Userdata');
            d = app.objects.STREAMobj.data{iy}.distance;
            A = app.objects.STREAMobj.data{iy}.IXgrid(d==0);
            linkobjname = app.objects.STREAMobj.names{iy};
            if ~isempty(A)
                % link distance to lowest point of streamobj
                itemname = app.objects.(classname).names{ix};
                x = app.objects.(classname).data{ix}.x;
                y = app.objects.(classname).data{ix}.y;
                iz = coord2ind(app.DEM,x,y);
                iz = iz(~isnan(iz));
                iz_fa = app.FA.Z(iz);
                iz = iz(iz_fa==max(iz_fa));
                [B,dist] = flowpathextract(app.FD,iz(1));
                [~,~,ib] = intersect(A,B,'stable');
                if ~isempty(ib)
                    offset = dist(ib);
                    fprintf(1,'Distance of %s linked to lowest point of %s.\n',itemname,linkobjname);
                else
                    fprintf(1,'Could not link %s to STREAMobj %s.\n',itemname,linkobjname);
                end
            end
        end
    end % 
    %------------------------------------------------------------
    
    %------------------------------------------------------------
    function mkswath(~,~,app)
        
        [classname,~,ix] = getobjectitem(app);
        if ~isempty(ix) && length(ix)==1
            
            hfig = gcf;
            busypointer(app,1)
            
            % resample dx
            h = findobj(hfig,'String','Resample x (m)');
            if get(h,'Value'); dx = str2double(get(get(h,'Userdata'),'String'));
            else dx = []; end
            
            % resample dy
            h = findobj(gcf,'String','Resample y (m)');
            if get(h,'Value'); dy = str2double(get(get(h,'Userdata'),'String'));
            else dy = []; end
            
            % keep nodes
            h = findobj(gcf,'String','Keep nodes');
            if get(h,'Value'); keepnodes = true; else keepnodes = false; end
            
            % smoothing
            h = findobj(gcf,'String','Smoothing (set below)');
            if get(h,'Value'); smooth = round(get(get(h,'Userdata'),'Value'));
            else smooth = 0; end
                        
            % swath width
            h = findobj(gcf,'String','Swath width (m)');
            g = get(h,'Userdata');
            width = str2double(get(g,'String'));
            
            % central gap
            h = findobj(gcf,'String','Central gap (m)');
            g = get(h,'Userdata');
            gap = max(0,str2double(get(g,'String')));
            
            % keep distance
            h = findobj(gcf,'String','Keep distance');
            if get(h,'Value'); keepdist = true; else keepdist = false; end
            
            % create SWATHobj
            if strcmp(classname,'STREAMobj');
                S = app.objects.(classname).data{ix}; %%%
                fprintf(1,'TOPOAPP: Swath mapping from STREAMobj only applied to trunk stream\n');
                S = trunk(S);
                SW = STREAMobj2SWATHobj(S,app.DEM,'dx',dx,'dy',dy,'width',...
                    width,'gap',gap,'smooth',smooth,'smoothlongest',true,...
                    'keepnodes',keepnodes,'keepdist',keepdist,'plot',false);
            else
                x = app.objects.(classname).data{ix}.x;
                y = app.objects.(classname).data{ix}.y;
                XY = [x,y];
                SW = SWATHobj(app.DEM,XY,'dx',dx,'dy',dy,'width',width,...
                    'gap',gap,'smooth',smooth,'smoothlongest',true,...
                    'keepnodes',keepnodes,'keepdist',keepdist,'plot',false);
            end
            
            % remove duplicates from mapping
            h = findobj(gcf,'String','Remove duplicates');
            if get(h,'Value'); SW = tidy(SW); end
            
            % Plot swath in map-view
            % Left/right
            h = findobj(gcf,'String','Left');
            if get(h,'Value'); left=true; else left=false; end
            h = findobj(gcf,'String','Right');
            if get(h,'Value'), right=true; else right=false; end
            % Outline
            h = findobj(gcf,'String','Show outline');
            if get(h,'Value'); outline=true; else outline=false; end
            % Points
            h = findobj(gcf,'String','Show points');
            if get(h,'Value'), points=true; else points=false; end
            
            axes(app.gui.hax), hold on
            SW.X = interp1(app.X,(1:app.DEM.size(2)),SW.X);
            SW.Y = interp1(app.Y,(1:app.DEM.size(1)),SW.Y);
            
            h = findobj('Tag','TemporaryObject');
            try delete(h); end
            h = plot(SW,'outline',outline,'points',points,'left',left,'right',right,'legend',false); hold off
            set(h,'Tag','TemporaryObject')
            figure(hfig)
            busypointer(app,0)
            
        else
            error('Topoapp: use only one item for swath mapping')
        end
    end %
    %------------------------------------------------------------
    function plotswath(hobject,~,app)
        
        hfig = gcf;
        [classname,~,ix] = getobjectitem(app);
        if isa(classname,'SWATHobj')
            SW = app.objects.(classname).data{ix};
            
        else
            h = findobj(app.gui.hax,'Tag','TemporaryObject');

            if isempty(h)
                warning('No swath data found. Generate swath with ''Map swath'' first');
                
            else
                
                distadjust = adjustdistance(classname,ix,app);
                if abs(distadjust)>0
                    % check if 'keep distance' activated
                    h = findobj('String','Keep distance');
                    if ~get(h,'Value')
                        set(0,'DefaultUIcontrolUnits','pixel');
                        button = questdlg('''Link to STREAMobj'' activated but ''keep distance'' not activated. Activate ''keep distance'' now (recommended)?','Warning','Yes','No','Yes');
                        set(0,'DefaultUIcontrolUnits','normalized');
                        if strcmp(button,'Yes'); set(h,'Value',1); end
                    end
                end
                

                switch get(hobject,'String')
                    case 'Plot elevation swath'
                        ylabelstr = 'Elevation (m)';
                        SW.zunit = 'Elevation (m)';
                    case 'Plot slope swath'
                        [SW] = mapswath(SW,app.G);
                        ylabelstr = 'Slope (deg)';
                        SW.zunit = 'Slope (deg)';
                end
                
                % Left/right
                h = findobj(gcf,'String','Left');
                if get(h,'Value'); left=true;
                else left=false; end
                h = findobj(gcf,'String','Right');
                if get(h,'Value'), right=true;
                else right=false; end

                figure
                plotdz(SW,'left',left,'right',right);
                ylabel(ylabelstr)
                updateaxis(gca)
                figure(hfig)
                busypointer(app,0)
                
            end
        end
    end % 
    %------------------------------------------------------------
    function [app] = addswath(hobject,~,app)
        if ~isempty(SW)
            axes(app.gui.hax), hold on
            h = plotobject(app,SW,app.objects.(class(SW)).color); hold off
            set(h,'Tag',class(SW))
            [app] = addobject(app,SW,'handle',h,'visible',true);
            set(h,'DisplayName',app.objects.(class(SW)).names{end})
        else
            warning('No SWATHobj available. Generate one by plotting first.')
        end
    end
    %------------------------------------------------------------

end %