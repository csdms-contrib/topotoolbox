function app = slopeareaplot(hObject,eventdata,app)
% SLOPEAREAPLOT provides access to TopoToolbox function 'slopearea'
%
% See also: topoapp/listfigure, topoapp/initclass, topoapp/addobject,
% slopearea, STREAMobj/modify
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2013

if strcmp(eventdata,'init') % initialize tool
    
    % Load button icon
    slopeareaicon = imread('slopeareaicon.png','png');
    
    % Set-up toolbar button
    app.gui.TB(end+1) = uipushtool('Parent',app.gui.hTB,...
        'Cdata',slopeareaicon,'TooltipString','Slope-area plot',...
        'ClickedCallback',{@slopeareaplot,app});
    
else % execute tool
    
    if isempty(app.S)
        warning('No STREAMobj found. Use FLOW routing button first')
    else
        % Initialize basic figure
        [hf,posYlast] = listfigure(app,{'STREAMobj'});
        % Add tool-specific options
        addoptions(hf,posYlast,app);
    end
end
    %------------------------------------------------------------
    function hfig = addoptions(hfig,posYlast,app)
        
        figure(hfig)
        %--------------------------------------------------------
        % positioning:
        dLs = 0.025; dLl = 0.2833;
        Lft1 = dLs; Lft2 = Lft1+(dLl+dLs+dLs); Lft3 = Lft2+(dLl+dLs+dLs);
        dLftedit1 = 0.12;  dLftedit2 = 0.2133;
        dLftpop1 = 0.155;
        Hpb = 0.06; Hedit = 0.045; Htxt = 0.045;
        Wpb = dLl; Wedit = 0.07; Wtxts = 0.1; Wtxtl = 0.17;
        Wpop1 = 0.13;
        %--------------------------------------------------------
        
        %% Parameter controls
        % Number of area bins
        posY = posYlast; ct = 1;
        h(ct) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,.15,Htxt],...
            'String','Area bins',...
            'TooltipString','Number of drainage area bins'); ct = ct+1;
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String','100',...
            'TooltipString','Number of drainage area bins');
        posY = posY - 0.05; ct = ct+1;
        % Minimum gradient
        h(ct) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,.15,Htxt],...
            'String','Min slope',...
            'TooltipString','Minimum slope value (m/m)'); ct = ct+1;
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String','1e-4',...
            'TooltipString','Minimum slope value (m/m)');
        posY = posY - 0.05; ct = ct+1;
        % Area bin locations
        h(ct) = uicontrol('Style','popupmenu',...
            'units','normalized','position',[Lft3+dLftpop1,posY,Wpop1,Hedit],...
            'String',{'median','mean','center'},...
            'Value',1); ct = ct+1;
         h(ct) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,.15,Htxt],...
            'String','Bin center',...
            'Userdata',h(ct-1),...
            'TooltipString','Determines the area values of the bin centers'); ct = ct+1;
        posY = posY - 0.05; 
        % Slope value calculation
        h(ct) = uicontrol('Style','popupmenu',...
            'units','normalized','position',[Lft3+dLftpop1,posY,Wpop1,Hedit],...
            'String',{'median','mean'},...
            'Value',1); ct = ct+1;
        h(ct) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,.15,Htxt],...
            'String','Slope',...
            'Userdata',h(ct-1),...
            'TooltipString','Determines how slope values are calculated in each area bin'); ct = ct+1;
        posY = posY - 0.05; 
        % Fit method
        h(ct) = uicontrol('Style','popupmenu',...
            'units','normalized','position',[Lft3+dLftpop1,posY,Wpop1,Hedit],...
            'String',{'least squares','least absolute deviations'},...
            'Value',1); ct = ct+1;
        h(ct) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,.15,Htxt],...
            'String','Fitting method',...
            'Userdata',h(ct-1),...
            'TooltipString','Fitting method of the slope-area data'); ct = ct+1;
        posY = posY - 0.05; 
        % 2-D density plot
        h(ct) = uicontrol('Style','checkbox',...
            'String','2-D density plot',...
            'units','normalized','Position',[Lft3,posY,Wpb,Hedit],...
            'Value',0);
        posY = posY - 0.05; ct = ct+1;
        
        % 2-D density bins
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit1,posY,Wedit,Hedit],...
            'String','100'); ct = ct+1;
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String','100'); ct = ct+1;
        h(ct) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,Wtxts,Htxt],...
            'TooltipString','Number of bins used for 2-D histogram',...
            'String','[n,m] bins','Userdata',h([ct-1,ct-2]));
        posY = posY - 0.05; ct = ct+1;
        
        % Strahler order
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit1,posY,Wedit,Hedit],...
            'String','min'); ct = ct+1;
        h(ct) = uicontrol('Style','edit',...
            'units','normalized','position',[Lft3+dLftedit2,posY,Wedit,Hedit],...
            'String','max'); ct = ct+1;
        h(ct) = uicontrol('Style','text',...
            'units','normalized','position',[Lft3,posY,Wtxts,Htxt],...
            'TooltipString','Select by channels by Strahler order',...
            'String','Strahler','Userdata',h([ct-1,ct-2])); ct = ct+1;

        %% Tool options
        h(ct) = uicontrol('Style','pushbutton',...
            'units','normalized','position',[Lft2,.275,Wpb,Hpb],...
            'String','Plot',...
            'Callback',{@mkslopeareaplot,app}); ct = ct+1;

        set(findobj('Style','text'),'Backgroundcolor',[.8 .8 .8])
        set(findobj('Style','checkbox'),'Backgroundcolor',[.8 .8 .8])
        %--------------------------------------------------------
        
    end %
    %------------------------------------------------------------
    function mkslopeareaplot(~,~,app)
        [classname,~,ix] = getobjectitem(app);
        if ~isempty(ix)
            for i = 1 : length(ix) % loop over objects
                s = app.objects.(classname).data{ix(i)};
                name = app.objects.(classname).names{ix(i)};
                [s] = limit2strahler(s);
                if isempty(s)
                    fprintf(1,'Current limits imposed by strahler order create empty STREAMobj./n');
                else
                    
                    h = findobj('Style','edit','TooltipString','Number of drainage area bins');
                    nbins = str2double(get(h,'String'));
                    if isnan(nbins); nbins = 100; end
                    set(h,'String',num2str(nbins));
                    
                    h = findobj('Style','edit','TooltipString','Minimum slope value (m/m)');
                    mingrad = str2double(get(h,'String'));
                    if isnan(mingrad); mingrad = 1e-4; end
                    set(h,'String',num2str(mingrad));
                    
                    h = get(findobj('String','Bin center'),'Userdata');
                    options = get(h,'String');
                    bincenter = options{get(h,'Value')};
                    
                    h = get(findobj('String','Slope'),'Userdata');
                    options = get(h,'String');
                    slopevalfct = options{get(h,'Value')};
                    
                    h = get(findobj('String','Fitting method'),'Userdata');
                    fitmethod = get(h,'Value');
                    if fitmethod==1
                        fitmethod = 'ls';
                    elseif fitmethod==2
                        fitmethod = 'lad';
                    end
                    
                    h = findobj('Style','checkbox','String','2-D density plot');
                    
                    if get(h,'Value')
                        h = get(findobj('String','[n,m] bins'),'Userdata');
                        n = str2double(get(h(1),'String'));
                        if isnan(n); n = 100; end
                        set(h(1),'String',num2str(n));
                        m = str2double(get(h(2),'String'));
                        if isnan(m); m = 100; end
                        set(h(2),'String',num2str(m));
                        
                        figure
                        SA = slopearea(s,app.DEM,app.FA,'areabins',nbins,'areabinlocs',bincenter,...
                            'gradaggfun',slopevalfct,'fitmethod',fitmethod,'hist2',true,'histbins',[n,m],'mingradient',mingrad);
                    else
                        figure
                        SA = slopearea(s,app.DEM,app.FA,'areabins',nbins,'areabinlocs',bincenter,...
                            'gradaggfun',slopevalfct,'fitmethod',fitmethod,'mingradient',mingrad);
                    end
                    
                    updateaxis(gca)
                    title(name)
                    str = sprintf('S = %1.2f A ^{%1.3f}',SA.ks,SA.theta);
                    xlims = get(gca,'Xlim'); ylims = get(gca,'Ylim');
                    text(xlims(1)*10,ylims(1)*10,str);
                    fprintf(1,'%s: S = %1.2f A ^{%1.3f}\n',name,SA.ks,SA.theta);
                end
            end
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


end



