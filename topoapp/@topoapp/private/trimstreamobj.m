function app = trimstreamobj(hObject,eventdata,app)
% TRIMSTREAMOBJ allows modifying existing STREAMobj in topoapp
%
% See also: topoapp/initclass, topoapp/addobject, STREAMobj,
% STEAMobj/modify
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013


if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    trimicon = imread('trimicon.png','png');
    
    % Set up toolbar button
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',trimicon,'TooltipString','Trim drainage network',...
        'ClickedCallback',{@trimstreamobj,app});
    
    % add class
    [app] = initclass(app,'STREAMobj','r');
    
else % execute tool
    
    if isempty(app.S)
        warning('No STREAMobj found. Use FLOW routing button first')
    else
        SObj = [];
        % Initialize basic figure
        [hf,posYlast] = listfigure(app,'STREAMobj');
        % Add tool-specific options
        addoptions(hf,posYlast,app);
    end
    
end %


    %------------------------------------------------------------
    function addoptions(hfig,posYlast,app)

        figure(hfig);
        %--------------------------------------------------------
        % positioning:
        dLs = 0.025; dLl = 0.2833;
        Lft1 = dLs; Lft2 = Lft1+(dLl+dLs+dLs); Lft3 = Lft2+(dLl+dLs+dLs);
        dLftedit1 = 0.12;  dLftedit2 = 0.2133;
        Hpb = 0.06; Hedit = 0.045; Htxt = 0.045;
        Wpb = dLl; Wedit = 0.07; Wtxts = 0.1; Wtxtl = 0.17;
        %--------------------------------------------------------
        
        %% Parameter controls
        posY = posYlast; ct = 1;
%         % Link distances to other streamobj
%         h(ct) = uicontrol('Style','checkbox',...
%             'units','normalized','position',[Lft3,posY,Wpb,Htxt],...
%             'String','Link distance to STREAMobj',...
%             'Callback',{@link2streamobj,app});
%         posY = posY - 0.05; ct = ct+1;
        % Strahler order
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit1,posY,Wedit,Hedit],...
            'String','min'); ct = ct+1;
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String','max'); ct = ct+1;
        h(ct) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,Wtxts,Htxt],...
            'String','Strahler','Userdata',h([ct-2,ct-1]));
            posY = posY - 0.05; ct = ct+1;
        % Minimum distance
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String','max'); ct = ct+1;
        h(ct) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,2*Wtxts,Htxt],...
            'String','Minimum distance','Userdata',h(ct-1));
            posY = posY - 0.05; ct = ct+1;
        % Tributary to
        %     'tributaryto' instance of STREAMobj
        %     returns the stream network that is tributary to a stream (network)        
        %--------------------------------------------------------
        
        %% Tool options
        % Apply modification
        h(ct) = uicontrol('Style','pushbutton',...
            'units','normalized','position',[Lft2,.875,Wpb,Hpb],...
            'String','Apply modification',...
            'Callback',{@applysettings,app}); ct = ct+1;
        % Get trunk
        h(ct) = uicontrol('Style','pushbutton',...
            'units','normalized','position',[Lft2,.8,Wpb,Hpb],...
            'String','Get trunk',...
            'Callback',{@gettrunk,app}); ct = ct+1;
        h(ct) = uicontrol('Style','pushbutton',...
            'units','normalized','position',[Lft2,.725,Wpb,Hpb],...
            'String','Add STREAMobj to topoapp',...
            'Callback',{@addSTREAMobj,app}); ct = ct+1;
        
        set(h,'Interruptible','off'); clear h;
        set(findobj('Style','text'),'Backgroundcolor',[.8 .8 .8])
        set(findobj('Style','checkbox'),'Backgroundcolor',[.8 .8 .8])
        %--------------------------------------------------------
        
    end
    %------------------------------------------------------------
    function [app] = applysettings(~,~,app)
        [classname,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            SObj = app.objects.(classname).data{ix};
            [SObj] = limit2strahler(SObj);
            [SObj] = mindist(SObj);
            if ~isempty(SObj.ix)
            h = findobj('Tag','TemporaryObject');
            try delete(h); end
            axes(app.gui.hax), hold on
            h = plotobject(app,SObj,'w'); hold off
            set(h,'Tag','TemporaryObject')
            else
                warning('Current settings produced empty STREAMobj.')
            end
        end
    end
    %------------------------------------------------------------
    function [app] = gettrunk(~,~,app)
        [classname,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            SObj = app.objects.(classname).data{ix};
            SObj = trunk(SObj);
            axes(app.gui.hax), hold on
            h = plotobject(app,SObj,app.objects.(class(SObj)).color);
            hold off
            set(h,'Tag',class(SObj))
            [app] = addobject(app,SObj,'handle',h,'visible',true);
            set(h,'DisplayName',app.objects.(class(SObj)).names{end})
            figure, plotdz(SObj,app.DEM);
        end
    end %
    %------------------------------------------------------------
    function [SObj] = limit2strahler(SObj)
        hsto = findobj('String','Strahler','Visible','on');
        if ~isempty(hsto)
            hst = get(hsto,'Userdata');
            stmin = str2double(get(hst(1),'String'));
            if isnan(stmin);
                set(hst(1),'String','min'); 
            else
                SObj = modify(SObj,'streamorder',['>=',num2str(stmin)]); 
            end
            stmax = str2double(get(hst(2),'String')); 
            if isnan(stmax);
                set(hst(2),'String','max');
            else
                SObj = modify(SObj,'streamorder',['<=',num2str(stmax)]); 
            end
        end
    end % 
    %------------------------------------------------------------
    function [SObj] = mindist(SObj)
        if ~isempty(SObj.ix)
            h = findobj(gcf,'String','Minimum distance');
            hmd = get(h,'Userdata');
            v = str2double(get(hmd,'String'));
            if isnan(v); 
                set(hmd,'String','max'); 
            else
                SObj = modify(SObj,'distance',v);
            end
        end
    end
    %------------------------------------------------------------
    function [app] = addSTREAMobj(hobject,~,app)
        if ~isempty(SObj)
            h = findobj('Tag','TemporaryObject');
            try delete(h); end
            axes(app.gui.hax), hold on
            h = plotobject(app,SObj,app.objects.(class(SObj)).color); hold off
            set(h,'Tag',class(SObj))
            [app] = addobject(app,SObj,'handle',h,'visible',true);
            set(h,'DisplayName',app.objects.(class(SObj)).names{end})
        else
            warning('No STREAMobj available. Generate one by applying settings first.')
        end
    end
    %------------------------------------------------------------
    
end